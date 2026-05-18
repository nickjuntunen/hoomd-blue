// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __UPDATER_MUVT_2D_H__
#define __UPDATER_MUVT_2D_H__

/*! \file UpdaterMuVT2D.h
 *  \brief Grand-canonical (muVT) updater corrected for 2D simulations.
 *
 *  \section thermodynamics Thermodynamic rationale
 *
 *  In three dimensions the grand-canonical acceptance criterion for particle
 *  insertion is
 *
 *      P_acc = min(1, z · V / ((N+1) · kT) · exp(−β ΔU))
 *
 *  where V = Lx Ly Lz is the box volume and z is the fugacity with units of
 *  [energy / volume].
 *
 *  In two dimensions the analogous criterion replaces V by the box area A:
 *
 *      P_acc = min(1, z₂D · A / ((N+1) · kT) · exp(−β ΔU))
 *
 *  where A = Lx Ly and z₂D has units of [energy / area].
 *
 *  HOOMD-blue uses 2D boxes with Lz = 0, so the raw BoxDim::getVolume() call
 *  in the base class returns zero, making all insertion/removal moves
 *  impossible.  UpdaterMuVT2D overrides getEffectiveVolume() to return the
 *  area A = Lx Ly via BoxDim::getVolume(true), restoring detailed balance.
 *
 *  The Gibbs-ensemble volume move (area move in 2D) is also handled correctly:
 *  both partitions use getEffectiveVolume() consistently, so the log-ratio
 *  acceptance term log(A_new / A) and the isotropic box rescaling L ∝ A^(1/2)
 *  (already implemented in the base class for 2D, line 1090–1097) remain
 *  self-consistent.
 *
 *  \section usage Usage
 *
 *  Use UpdaterMuVT2D (Python: hpmc.update.MuVT2D) in place of UpdaterMuVT
 *  whenever the simulation box has Lz = 0 (ndimensions == 2).  The fugacity
 *  parameter must have units of [energy · area⁻¹] instead of
 *  [energy · volume⁻¹].
 */

#include "UpdaterMuVT.h"

namespace hoomd
    {
namespace hpmc
    {

/*!
 * Grand-canonical (muVT) updater for 2D simulations.
 *
 * Inherits all logic from UpdaterMuVT<Shape> and overrides only
 * getEffectiveVolume() so that the acceptance criterion for insertion,
 * removal, and Gibbs-ensemble area moves uses the box area A = Lx · Ly
 * rather than the zero 3D volume Lx · Ly · Lz.
 */
template<class Shape> class UpdaterMuVT2D : public UpdaterMuVT<Shape>
    {
    public:
    //! Inherit all constructors from UpdaterMuVT<Shape> unchanged.
    using UpdaterMuVT<Shape>::UpdaterMuVT;

    protected:
    //! Return the box area (Lx · Ly) for use in 2D acceptance criteria.
    //! BoxDim::getVolume(true) returns m_L.x * m_L.y regardless of Lz.
    Scalar getEffectiveVolume(const BoxDim& box) const override
        {
        return box.getVolume(true);
        }
    };

namespace detail
    {

//! Export UpdaterMuVT2D<Shape> to Python.
/*!
 *  \param m    The pybind11 module to export into (the _hpmc extension module).
 *  \param name Python-visible class name, e.g. "UpdaterMuVT2DConvexPolygon".
 *  \tparam Shape  HPMC shape type (ShapeConvexPolygon, ShapeSphere, …).
 *
 *  The exported class inherits all Python-accessible properties from
 *  UpdaterMuVT<Shape> (fugacity, max_volume_rescale, transfer_types, …) and
 *  is constructible with the same arguments.
 */
template<class Shape>
void export_UpdaterMuVT2D(pybind11::module& m, const std::string& name)
    {
    pybind11::class_<UpdaterMuVT2D<Shape>,
                     UpdaterMuVT<Shape>,
                     std::shared_ptr<UpdaterMuVT2D<Shape>>>(m, name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            std::shared_ptr<Trigger>,
                            std::shared_ptr<IntegratorHPMCMono<Shape>>,
                            unsigned int>());
    }

    } // namespace detail
    } // namespace hpmc
    } // namespace hoomd

#endif // __UPDATER_MUVT_2D_H__
