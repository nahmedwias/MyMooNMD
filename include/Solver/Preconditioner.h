#ifndef __PRECONDITIONER_H__
#define __PRECONDITIONER_H__

/** @brief an abstract base class to describe preconditioners
 * 
 */
template <class Vector>
class Preconditioner
{
  public:
    /// @brief use this object as a preconditioner
    virtual void apply(const Vector & z, Vector & r) const = 0;
    
    /** @brief Method to use this preconditioner in FGMRES.
     *
     * So far i and j get ignored and simply solve(z,r) is called.
     *
     * @param i Number of current iteration since last restart in FGMRES.
     * @param j Number of current iteration in FGMRES.
     * @param z The right hand side of the preconditioning.
     * @param r The obtained vector.
     */
    void apply(unsigned int i, unsigned int j, const Vector &z, Vector &r) const
    { apply(z, r); }
};

template <class Vector>
class NoPreconditioner : public Preconditioner<Vector>
{
  public:
    
    void apply(const Vector & z, Vector & r) const override
    { r = z; }
    
    /** @brief Method to use this preconditioner in FGMRES.
     *
     * So far i and j get ignored and simply solve(z,r) is called.
     *
     * @param i Number of current iteration since last restart in FGMRES.
     * @param j Number of current iteration in FGMRES.
     * @param z The right hand side of the preconditioning.
     * @param r The obtained vector.
     */
    void apply(unsigned int i, unsigned int j, const Vector &z, Vector &r) const
    { apply(z, r); }
};

#endif // __PRECONDITIONER_H__