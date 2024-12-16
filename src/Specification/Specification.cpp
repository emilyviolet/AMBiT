#include <algorithm>
#include <fstream>
#include <istream>
#include <string>
#include <string_view>

#include "Specification.h"
#include "parapara/parapara.h"

namespace P = parapara;

namespace Ambit
{

template <typename Record>
P::hopefully<void> custom_import_ini(Record& rec, const P::specification_set<Record>& specs, std::istream& in);

inline auto positive = P::greater_than(0, "must be positive");

P::specification<GlobalSpecification> global_specifications[] = {
    {"Lattice/NumPoints",     &GlobalSpecification::lattice_num_points,    P::nonzero()},
    {"Lattice/StartPoint",    &GlobalSpecification::lattice_start_point,   positive},
    {"Lattice/EndPoint",      &GlobalSpecification::lattice_end_point,     positive},
    {"Lattice/--exp-lattice", &GlobalSpecification::lattice_exponential},
    {"Lattice/H",             &GlobalSpecification::lattice_H,             positive}
};

P::specification_set<GlobalSpecification> global_specifications_dict(global_specifications, P::keys_lc_nows);

LatticeConfig GlobalSpecification::getLatticeConfig() const {
    if (lattice_exponential) {
        LatticeExpConfig config;

        config.H = lattice_H;
        if (lattice_num_points>0) config.num_points = lattice_num_points;
        if (lattice_start_point>0) config.start_point = lattice_start_point;
        return config;
    }
    else {
        LatticeHybridConfig config;

        if (lattice_num_points>0) config.num_points = lattice_num_points;
        if (lattice_start_point>0) config.start_point = lattice_start_point;
        if (lattice_end_point>0) config.end_point = lattice_end_point;
        return config;
    }
}

std::string importSpecificationFile(GlobalSpecification& gs, const std::string& fileName) {
    std::ifstream in(fileName);
    if (!in) return "unable to open input file '"+fileName+"'";

    auto h = custom_import_ini(gs, global_specifications_dict, in);
    if (!h) return P::explain(h.error(), true);

    return "";
}

std::string importSpecificationKV(GlobalSpecification& gs, const std::string& assignment) {
    auto h = P::import_k_eq_v(gs, global_specifications_dict, assignment);
    if (!h) return P::explain(h.error(), true);
    else return "";
}

// Perform global validation of specifications. Return non-empty error message on failure.
std::string validateSpecification(const GlobalSpecification& gs) {
    // Check for consistent Lattice settings:

    if (gs.lattice_exponential && gs.lattice_end_point>0)
        return "Lattice/EndPoint is ignored if Lattice/--exp-lattice is set";

    return "";
}

// Custom INI import implementation:

using P::ini_record;
using P::ini_record_kind;

ini_record custom_ini_parser(std::string_view v) {
    using token = ini_record::token;
    using size_type = std::string_view::size_type;

    constexpr size_type npos = std::string_view::npos;
    constexpr std::string_view ws{" \t\f\v\r\n"};

    auto comment_at = [](std::string_view v, size_type p) { return v.substr(p, 2)=="//"; };
    size_type b = v.find_first_not_of(ws);

    // empty or comment?
    if (b==npos || comment_at(v, b)) return ini_record{ini_record_kind::empty};

    // section heading?
    if (v[b]=='[') {
        size_type e = v.find(']');

        // check for malformed heading
        if (e==npos) return {ini_record_kind::syntax_error, token("", b+1)};

        if (e+1<v.length()) {
            auto epilogue = v.find_first_not_of(ws, e+1);
            if (epilogue!=npos && !comment_at(v, epilogue)) {
                return {ini_record_kind::syntax_error, token("", epilogue)};
            }
        }

        b = v.find_first_not_of(ws, b+1);
        e = v.find_last_not_of(ws, e-1);

        return {ini_record_kind::section, token(v.substr(b, e>=b? e+1-b: 0), b+1)};
    }

    // expect key first, followed by ws and eol, =, or //.
    size_type j = std::min(v.find('=', b), v.find("//", b));
    token key_token{v.substr(b, j==b? 0: v.find_last_not_of(ws, j-1)+1-b), b+1};

    // key without value?
    if (j==npos || v[j]!='=') {
        return {ini_record_kind::key, key_token};
    }

    // skip to text after =, look for value
    size_type eq = j;
    size_type value_cindex = eq;

    if (j<v.length()) {
        j = v.find_first_not_of(ws, j+1);
        if (j!=npos && !comment_at(v, j)) {
            value_cindex = j+1;

            // if value is not quoted, take text up to eol or first eol, discarding trailing ws
            if (v[j]!='\'') {
                size_type end = v.find("//", j);
                if (end!=npos) --end;

                return {ini_record_kind::key_value, key_token,
                        token{v.substr(j, v.find_last_not_of(ws, end)-j+1), value_cindex}};
            }
            else {
                // quoted value; take until next unescaped '
                std::string value;
                size_type epilogue = npos;
                bool esc = false;

                for (size_type i = j+1; i<v.length(); ++i) {
                    if (esc) {
                        value += v[i];
                        esc = false;
                    }
                    else if (v[i]=='\'') {
                        epilogue = i+1;
                        break;
                    }
                    else if (v[i]=='\\') {
                        esc = true;
                    }
                    else {
                        value += v[i];
                    }
                }

                // unterminated quoted value or escaped eol?
                if (epilogue==npos || esc) {
                    return {ini_record_kind::syntax_error, token("", j)};
                }

                // extra stuff following value that is not a comment?
                if (epilogue<v.length()) {
                    epilogue = v.find_first_not_of(ws, epilogue);
                    if (epilogue!=npos && !comment_at(v, epilogue)) {
                        return {ini_record_kind::syntax_error, token("", epilogue)};
                    }
                }

                return {ini_record_kind::key_value, key_token, token{value, value_cindex}};
            }
        }
    }
    // key with empty value
    return {ini_record_kind::key_value, key_token, token{"", eq}};
}

// Use a custom line-by-line ini importer to handle relative section headings
template <typename Record>
P::hopefully<void> custom_import_ini(Record& rec, const P::specification_set<Record>& specs, std::istream& in)
{
    constexpr auto npos = std::string_view::npos;

    P::ini_style_importer importer(custom_ini_parser, in);
    while (importer) {
        std::string prev_sec{importer.section()};
        auto h = importer.run_one(rec, specs, P::default_reader(), "/");

        // Ignore unrecognized keys until specification coverage is complete
        if (!h) {
             if (h.error().error!=P::failure::unrecognized_key) return P::unexpected(h.error());
             else continue; 
        }

        if (h.value() != P::ini_record_kind::section) continue;

        std::string_view new_sec = importer.section();
        if (new_sec.substr(0, 2)=="./") {
            std::string section{prev_sec};
            section += new_sec.substr(1);
            importer.section(section);
        }
        else if (new_sec.substr(0, 3)=="../") {
            std::string_view prefix(prev_sec);
            do {
                new_sec.remove_prefix(3);
                if (auto tail = prefix.rfind('/'); tail != npos) {
                    prefix = prefix.substr(0, tail);
                }
            } while (new_sec.substr(0, 3)=="../");

            std::string section{prefix};
            section += "/";
            section += new_sec;
            importer.section(section);
        }
    }
    return {};
}

} // namespace Ambit
