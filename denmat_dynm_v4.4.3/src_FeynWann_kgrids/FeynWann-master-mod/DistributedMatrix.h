/*-------------------------------------------------------------------
Copyright 2018 Ravishankar Sundararaman

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#ifndef FEYNWANN_DISTRIBUTEDMATRIX_H
#define FEYNWANN_DISTRIBUTEDMATRIX_H

#include <core/matrix.h>

//! Distributed matrix elements / Hamiltonians for FeynWann
class DistributedMatrix
{
public:
	//Input parameters:
	const MPIUtil* mpiUtil; //!< MPI instance / communicator over whcih this is parallelized
	int nElemsTot; //!< total number of matrix elements (per cell/k) [if packed, only counts packed]
	const std::vector<vector3<int>>& cellMap; //!< cell map in real space
	const vector3<int>& kfold; //!< k-point folding
	bool squared; //!< whether cell map / kpoints are squared (eg. e-ph matrix elements)
	
	//Parameters set in the constructor:
	int nElems; //!< local number of matrix elements (per cell/k) [if packed, only counts packed]
	int iElemStart; //!< starting element on current process
	int nCellsTot; //!< size of cell map or its square, depending on squared
	int kfoldProd; //!< prod(kfold)
	int nkTot; //!< prod(kfold) or its square, depending on squared
	int nk; //!< number of k-points (or pairs, if squared) on current process
	int ikStart; //!< starting k-point (or pair, if squared) on current process
	std::vector<int> ikStartProc; //!< ikStart for each process in MPI group
	std::vector<int> iElemStartProc; //!< iElemStart for each process in MPI group
	
	//! Initialize from file, containing complex or real elements as specified by realOnly
	//! If mpiInterGroup is non-null, then read only in one group and broadcast to rest (require regular process grid)
	//! If cellWeights is non-null, then read in only unique cells and use cellWeights in interpolation (required in squared mode)
	//! (remaining parameters are as specified in the class)
	DistributedMatrix(string fname, bool realOnly, const MPIUtil* mpiUtil, int nElemsTot,
		const std::vector<vector3<int>>& cellMap, const vector3<int>& kfold, bool squared,
		const std::shared_ptr<MPIUtil> mpiInterGroup=0, const std::vector<matrix>* cellWeights=0);
	~DistributedMatrix();
	
	void transform(vector3<> k0); //!< prepare results for k-point mesh offset by k0 (squared=false only)
	void transform(vector3<> k01, vector3<> k02); //!< prepare results for k-point mesh offsets k01 and k02 (squared=true only)
	const complex* getResult(int ik) const; //!< get pointer to result for k-point (or pair) index ik
	
	//Alternatve interface for single k-points (insted of offsets):
	void compute(vector3<> k); //!< prepare results for single point k (squared=false only)
	void compute(vector3<> k1, vector3<> k2); //!< prepare results for single pair k1,k2 (squared=true only)
	
private:
	ManagedArray<complex> mat; //!< input matrix elements
	ManagedArray<complex> buf; //!< buffer in which transformations happen and result is produced
	std::shared_ptr<struct PlanSet> planSet; //!< opaque pointer to required set of FFT plans
	std::vector<int> cellIndex; //!< index of cell in nkTot array
	//Cells and weights by unique indices:
	struct Cell
	{	vector3<int> iR;
		std::vector<double> weight; //nBands x nBands in non-squared, and nAtoms (outer) by nBands (inner) in squared mode
		complex phase01, phase02; //temporary variables used in transform
	};
	std::vector<std::vector<Cell>> uniqueCells;
	int nAtoms, nBands;
	void collectHead(); //!< collect buf results on mpiUtil head for compute (single k-point versions of transform)
};

//Calculate flat index given 3D coordinates and sample counts
inline int calculateIndex(const vector3<int>& iv, const vector3<int>& S)
{	int i = 0;
	for(int iDir=0; iDir<3; iDir++)
	{	if(iDir) i *= S[iDir];
		i += positiveRemainder(iv[iDir], S[iDir]);
	}
	return i;
}

#endif //FEYNWANN_DISTRIBUTEDMATRIX_H
