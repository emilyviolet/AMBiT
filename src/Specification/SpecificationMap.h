#ifndef SPECIFICATION_MAP_H
#define SPECIFICATION_MAP_H

#include <any>
#include <string_view>
#include <type_traits>
#include "parapara/parapara.h"
struct record_field_proxy {
    explicit record_field_proxy(parapara::any_ptr p, std::function<void (std::any)> setter = {}):
        p_(std::move(p)), setter_(std::move(setter)) {}

    // Restrict methods to rvalues so that the proxy object cannot easily be used
    // outside the expression in which it is instantiated.

    template <typename X> const X& as() const && { return deref<X>(); }
    template <typename X> operator X() const && { return deref<X>(); }

    template <typename X> record_field_proxy&& operator=(X&& x) && {
        if (!setter_) {
            using Y = std::remove_const_t<std::remove_reference_t<X>>;
            if (Y* yp = p_.as<Y*>()) *yp = std::forward<X>(x);
            else throw std::bad_any_cast();
        }
        else {
            // Use supplied setter (so that we can hook in specification validators):
            setter_(std::forward<X>(x));
        }
        return std::move(*this);
    }

private:
    parapara::any_ptr p_;
    std::function<void (std::any)> setter_;

    template <typename X> const X& deref() const {
        if (X* xp = p_.as<X*>()) return *xp;
        if (const X* xp = p_.as<const X*>()) return *xp;
        throw std::bad_any_cast();
    }
};

// Keyed record view: hold on to a copy of a specification_map
// and a (non-owning) pointer to a corresponding record instance; provide
// interface to access fields of the record via the keys in
// the specification_map.

// Record-parameterized version:

template <typename Record = void>
struct keyed_record_view {
    explicit keyed_record_view(const parapara::specification_map<Record>& smap):
        smap_(smap) {}

    explicit keyed_record_view(const parapara::specification_map<Record>& smap, Record& record):
        smap_(smap), rptr_(&record) {}

    keyed_record_view& rebind(Record& record) {
        rptr_ = &record;
        return *this;
    }

    record_field_proxy operator[](std::string_view key) {
        // Throw on absent key.
        const parapara::specification<Record>& spec = smap_.at(key);

        if (!rptr_) {
            return record_field_proxy(parapara::any_ptr{});
        }
        else {
            // (Handling optional members currently a bit clunky; can be improved in parapara.)
            parapara::any_ptr p = spec.retrieve(*rptr_).value_or(nullptr); // empty any_ptr if empty optional field.
            return record_field_proxy(p,
                [spec, r = rptr_](std::any a) {
                    if (auto h = spec.assign(*r, std::move(a))) return;
                    else {
                        // assign() gives an 'internal_error' on type mismatch because this use-case wasn't anticipated.
                        // We can extend parapara to add a new error enum or reuse 'unsupported_type' as we do here.
                        parapara::failure fail(std::move(h).error());
                        if (fail.error==parapara::failure::internal_error) fail.error = parapara::failure::unsupported_type;
                        throw parapara::parapara_error(explain(fail));
                    }
                });
        }
    }

private:
    parapara::specification_map<Record> smap_;
    Record* rptr_ = nullptr;
};


// Generic version TODO

template <> struct keyed_record_view<void> {
};
#endif // SPECIFICATION_MAP_H
