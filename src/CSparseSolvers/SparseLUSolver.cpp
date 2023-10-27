/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define CSPARSESOLVERS_SPARSELUSOLVER_CPP
#include <CSparseSolvers/config.h>
#include <CSparseSolvers/SparseLUSolver.inl>
#include <sofa/core/ObjectFactory.h>

namespace csparsesolvers
{

using namespace sofa::linearalgebra;

// dont try to compile if floating point is float
#ifdef SOFA_FLOAT
SOFA_PRAGMA_WARNING("SparseLUSolver does not support float as scalar.")
#else  // SOFA_DOUBLE
int SparseLUSolverClass = sofa::core::RegisterObject("Direct linear solver based on Sparse LU factorization, implemented with the CSPARSE library")
        .add< SparseLUSolver< CompressedRowSparseMatrix<SReal>, FullVector<SReal> > >()
        .add< SparseLUSolver< CompressedRowSparseMatrix<sofa::type::Mat<3,3,SReal> >,FullVector<SReal> > >()
        ;

template class SOFA_CSPARSESOLVERS_API SparseLUSolver< CompressedRowSparseMatrix<SReal>, FullVector<SReal> >;
template class SOFA_CSPARSESOLVERS_API SparseLUSolver< CompressedRowSparseMatrix<sofa::type::Mat<3,3,SReal> >, FullVector<SReal> >;
#endif  // SOFA_FLOAT

} // namespace sofa::component::linearsolver::direct

