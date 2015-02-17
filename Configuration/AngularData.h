#ifndef ANGULAR_DATA_H
#define ANGULAR_DATA_H

#include "HartreeFock/OrbitalInfo.h"
#include "Projection.h"
#include "Symmetry.h"
#include <list>
#include <unordered_map>
#include <memory>
#include <boost/functional/hash.hpp>
#include <boost/iterator/iterator_facade.hpp>

class RelativisticConfiguration;
class AngularDataLibrary;

/** Store projections and configuration state functions (CSF) corresponding to
    the angular part of a RelativisticConfiguration.
    Access is only provided by const iterators.
 */
class AngularData
{
protected:
    friend class AngularDataLibrary;
    AngularData(int two_m);

public:
    AngularData(const RelativisticConfiguration& config, int two_m);    //!< Generate projections but not CSFs.
    AngularData(const RelativisticConfiguration& config, int two_m, int two_j); //!< Generate projections and CSFs.
    ~AngularData();

    typedef std::vector< std::pair<int, int> > ConfigKeyType;   // pair(kappa, number of particles) for all orbitals in RelativisticConfiguration
    typedef std::list< std::vector<int> >::const_iterator const_projection_iterator;
    typedef const double* const_CSF_iterator;

    /** Bidirectional list iterator over list of "projections", i.e. list of std::vector<int> */
    const_projection_iterator projection_begin() const { return projections.begin(); }
    const_projection_iterator projection_end() const { return projections.end(); }
    unsigned int projection_size() const { return projections.size(); }

    /** Return whether CSFs have been calculated (or read in). */
    bool CSFs_calculated() const { return have_CSFs; }

    /** Random access iterator over CSFs corresponding to projection i.
        PRE: 0 <= i <= projection_size()
     */
    const_CSF_iterator CSF_begin(int i) const { return CSFs + i * num_CSFs; }
    const_CSF_iterator CSF_end(int i) const { return CSFs + (i+1) * num_CSFs; }

    int GetTwoM() const { return two_m; }
    int GetTwoJ() const { return two_j; }
    int NumCSFs() const { return num_CSFs; }
    const double* GetCSFs() const { return CSFs; }

    /** Generate CSFs by diagonalising projections over J^2.
        Return number of CSFs generated with correct two_j.
     */
    int GenerateCSFs(const RelativisticConfiguration& config, int two_j);
    int GenerateCSFs(const ConfigKeyType& key, int two_j);

protected:
    int GenerateProjections(const RelativisticConfiguration& config, int two_m);
    static bool ProjectionCompare(const std::vector<int>& first, const std::vector<int>& second);

    /** List of "projections": in this context, vectors of two_Ms. */
    std::list< std::vector<int> > projections;
    int two_m;

    /** CSF coefficients for a given J. Usually one requires all coefficients for a given projection,
        so the projection comes first: CSFs[proj * N + csf] where N is NumCSFs().
     */
    bool have_CSFs;
    double* CSFs;
    int num_CSFs;
    int two_j;

protected:
    /** J^2 = Sum_i j_i^2 + 2 Sum_(i<j) j_i . j_j
        One-body operator = j(j+1)
        Two-body operator = 2 (jz_i.jz_j) + (j+_i.j-_j + j-_i.j+_j)
     */
    class JSquaredOperator
    {
    public:
        double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const
        {
            if(e1 == e2)
            {
                double value = e1.J() * (e1.J() + 1);
                // Even the one-body part of J^2 is not a real one body operator
                if(e1.IsHole())
                    return -value;
                else
                    return value;
            }
            else
                return 0.0;
        }

        double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const
        {
            if(e1.Kappa() == e3.Kappa() && e1.PQN() == e3.PQN()
               && e2.Kappa() == e4.Kappa() && e2.PQN() == e4.PQN())
            {
                // 2 (jz_i.jz_j)
                if(e1.TwoM() == e3.TwoM() && e2.TwoM() == e4.TwoM())
                    return 2. * e1.M() * e2.M();

                // (j+_i.j-_j + j-_i.j+_j)
                if((abs(e1.TwoM() - e3.TwoM()) == 2) && (abs(e2.TwoM() - e4.TwoM()) == 2)
                    && (e1.TwoM() + e2.TwoM() == e3.TwoM() + e4.TwoM()))
                {
                    return (sqrt(e1.J() * (e1.J() + 1.) - e1.M() * e3.M()) *
                            sqrt(e2.J() * (e2.J() + 1.) - e2.M() * e4.M()));
                }
            }

            return 0.;
        }
    } J_squared_operator;
};

typedef std::shared_ptr<AngularData> pAngularData;
typedef std::shared_ptr<const AngularData> pAngularDataConst;

/** Collection of AngularData elements, indexed by RelativisticConfiguration.
    The collection is stored on disk in the directory specified by lib_directory, with filename
        <particle_number>.<two_j>.<parity>.two_m.angular
    As usually specified in the Makefile, the directory is
        AMBiT/AngularData/
 */
class AngularDataLibrary
{
public:
    /** Initialise with number of electrons, symmetry, and directory of stored libraries.
        If lib_directory is not specified then Read() and Write() are disabled, so all CSFs must be recalculated.
        If lib_directory does not exist already, then this is an error so that the user can check the path specification.
     */
    AngularDataLibrary(int electron_number, const Symmetry& sym, int two_m, std::string lib_directory = "");
    ~AngularDataLibrary() {}

    /** Retrieve or create an AngularData object for the given configuration, two_m, and two_j.
        If there are no projections possible with our two_m, it returns a null pointer.
        PRE: config is sorted correctly.
     */
    pAngularData operator[](const RelativisticConfiguration& config);

    void GenerateCSFs();

    void Read();
    void Write() const;

protected:
    typedef AngularData::ConfigKeyType KeyType;

    int electron_number, two_m, two_j;
    KeyType GenerateKey(const RelativisticConfiguration& config) const;

    std::string filename;
    std::unordered_map<KeyType, pAngularData, boost::hash<KeyType> > library;
};

typedef std::shared_ptr<AngularDataLibrary> pAngularDataLibrary;

#endif
