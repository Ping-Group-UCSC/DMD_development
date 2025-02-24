#include "lindblad/LindbladFile.h"
#include "InputMap.h"
#include <core/Units.h>
#include <core/Minimize.h>
#include <core/Random.h>
#include <deque>

struct TripletMatrix : public LinearSolvable<matrix>
{
	//Triplet format:
	const int nRows, nCols;
	const bool symmetric; //whether the matrix is symmetric (implicit Mji for every Mij added)
	struct Entry
	{	int i,j; double Mij; 
		Entry(int i, int j, double Mij) : i(i), j(j), Mij(Mij) {}
	};
	std::vector<Entry> entries; //MPI divided; each process has a subset
	std::vector<bool> included; //whether each row/col index is included in matrix (all true by default)
	diagMatrix K; //preconditioner (identity by default)
	
	TripletMatrix(int nRows, int nCols, int nNZestimate=0, bool symmetric=false)
	: nRows(nRows), nCols(nCols), symmetric(symmetric), included(nRows, true), K(nRows, 1.)
	{	if(nNZestimate) entries.reserve(nNZestimate / mpiWorld->nProcesses());
	}
	
	//Matrix multiply:
	matrix operator*(const matrix& v) const
	{	static StopWatch watch("TripletMatrix::operator*"); watch.start();
		matrix out = zeroes(nRows, v.nCols());
		complex* outData = out.data(); const complex* vData = v.data();
		for(const Entry& entry: entries)
			for(int col=0; col<v.nCols(); col++)
			{	outData[out.index(entry.i, col)] += entry.Mij * vData[v.index(entry.j, col)];
				if(symmetric)
					outData[out.index(entry.j, col)] += entry.Mij * vData[v.index(entry.i, col)];
			}
		mpiWorld->allReduceData(out, MPIUtil::ReduceSum);
		watch.stop();
		return out;
	}
	inline matrix hessian(const matrix& v) const { return (*this) * v; }
	
	//Calculate Sum_i Mij:
	diagMatrix DiagSumRows() const
	{	diagMatrix out(nCols, 0.);
		for(const Entry& entry: entries)
		{	out[entry.j] += entry.Mij;
			if(symmetric) out[entry.i] += entry.Mij;
		}
		mpiWorld->allReduceData(out, MPIUtil::ReduceSum);
		return out;
	}
	
	double sync(double x) const { mpiWorld->bcast(x); return x; }
	
	//Solve A x = b iteratively using CG
	void applyInverse(const matrix& b, matrix& x, double relTol=1e-9)
	{	MinimizeParams mp;
		mp.nIterations = 1000;
		mp.fpLog = nullLog;
		for(int col=0; col<b.nCols(); col++)
		{	state = project(x(0,nRows, col,col+1));
			matrix bCol = project(b(0,nRows, col,col+1));
			mp.knormThreshold = relTol * sqrt(dot(bCol, precondition(bCol)));
			logPrintf("\tCG solve: "); logFlush();
			int nIter = solve(bCol, mp);
			logPrintf("%s after %d iterations.\n", (nIter<mp.nIterations ? "converged" : "did not converge"), nIter); logFlush();
			x.set(0,nRows, col,col+1, state);
		}
	}
	
	inline matrix project(const matrix& v) const
	{	matrix out(v);
		complex* outData = out.data();
		for(int col=0; col<out.nCols(); col++)
		{	complex mean = 0.;
			int nContrib = 0;
			for(int row=0; row<nRows; row++)
				if(included[row])
				{	mean += outData[out.index(row,col)];
					nContrib++;
				}
			mean *= (1./nContrib);
			for(int row=0; row<nRows; row++)
			{	if(included[row])
					outData[out.index(row,col)] -= mean;
				else 
					outData[out.index(row,col)] = 0.;
			}
		}
		return out;
	}
	
	matrix precondition(const matrix& v) const
	{	return project(K * project(v));
	}
};

matrix3<> Matrix3(const matrix& m)
{	assert(m.nRows()==3);
	assert(m.nCols()==3);
	matrix3<> result;
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			result(i,j) = m(i,j).real();
	return result;
}

int main(int argc, char** argv)
{	
	InitParams ip = FeynWann::initialize(argc, argv, "Linearized-Boltzmann calculation of resistivity");

	//Read input file:
	InputMap inputMap(ip.inputFilename);
	const double Tmin = inputMap.get("Tmin") * Kelvin; //temperature; start of range
	const double Tmax = inputMap.get("Tmax", Tmin/Kelvin) * Kelvin; assert(Tmax>=Tmin); //temperature; end of range (defaults to Tmin)
	const size_t Tcount = inputMap.get("Tcount", 1); assert(Tcount>0); //number of temperatures
	const double dmu = inputMap.get("dmu", 0.) * eV; //must be within [dmuMin,dmuMax] specified while generating inFile
	const string inFile = inputMap.has("inFile") ? inputMap.getString("inFile") : "ldbd.dat"; //input file name
	
	logPrintf("\nInputs after conversion to atomic units:\n");
	logPrintf("Tmin = %lg\n", Tmin);
	logPrintf("Tmax = %lg\n", Tmax);
	logPrintf("Tcount = %lu\n", Tcount);
	logPrintf("dmu = %lg\n", dmu);
	logPrintf("inFile = %s\n", inFile.c_str());
	logPrintf("\n");
	
	//Read inFile generated by lindbladInit:
	logPrintf("Reading '%s': ", inFile.c_str()); logFlush();
	//--- Read header and check parameters:
	MPIUtil::File fp;
	mpiWorld->fopenRead(fp, inFile.c_str());
	LindbladFile::Header h; h.read(fp, mpiWorld);
	if(dmu<h.dmuMin or dmu>h.dmuMax)
	{	FeynWann::finalize();
		die("dmu = %lg eV is out of range [ %lg , %lg ] eV specified in lindbladInit.\n", dmu/eV, h.dmuMin/eV, h.dmuMax/eV);
	}
	if(Tmax > h.Tmax)
	{	FeynWann::finalize();
		die("Tmax = %lg K is larger than Tmax = %lg K specified in lindbladInit.\n", Tmax/Kelvin, h.Tmax/Kelvin);
	}
	double Omega = fabs(det(h.R)); //unit cell volume
	if(not h.ePhEnabled)
	{	FeynWann::finalize();
		die("resistivityBTE requires lindbladInit to be run with e-ph enabled.\n");
	}
	std::vector<size_t> byteOffsets(h.nk);
	mpiWorld->freadData(byteOffsets, fp);
	//--- Read energies and matrix elements divided by k-points:
	TaskDivision kDivision(h.nk, mpiWorld);
	size_t ikStart, ikStop;
	kDivision.myRange(ikStart, ikStop);
	size_t nkMine = ikStop-ikStart;
	size_t ikInterval = std::max(1, int(round(nkMine/50.))); //interval for reporting progress
	std::vector<LindbladFile::Kpoint> kpoint(nkMine);
	std::vector<int> nInnerAll(h.nk); //number of active bands at each k
	mpiWorld->fseek(fp, byteOffsets[ikStart], SEEK_SET);
	for(size_t ik=ikStart; ik<ikStop; ik++)
	{	LindbladFile::Kpoint& kp = kpoint[ik-ikStart];
		kp.read(fp, mpiWorld, h);
		nInnerAll[ik] = kp.nInner;
		//Print progress:
		if((ik-ikStart+1)%ikInterval==0) { logPrintf("%d%% ", int(round((ik-ikStart+1)*100./nkMine))); logFlush(); }
	}
	mpiWorld->fclose(fp);
	logPrintf("done.\n"); logFlush();
	logPrintf("%lu active k-points parallelized over %d processes.\n", h.nk, mpiWorld->nProcesses());
	
	if(ip.dryRun)
	{	logPrintf("Dry run successful: commands are valid and initialization succeeded.\n");
		FeynWann::finalize();
		return 0;
	}
	logPrintf("\n");

	//Determine total number of bands by k on all processes:
	for(int jProc=0; jProc<mpiWorld->nProcesses(); jProc++)
		mpiWorld->bcast(nInnerAll.data()+kDivision.start(jProc),
			kDivision.stop(jProc)-kDivision.start(jProc), jProc);
	std::vector<size_t> stateOffset(h.nk+1); //global index of first active band within each k
	for(size_t ik=0; ik<h.nk; ik++)
		stateOffset[ik+1] = stateOffset[ik] + nInnerAll[ik];
	size_t stateOffsetMine = stateOffset[ikStart]; //offset of states stored in current process
	size_t nStatesMine = stateOffset[ikStop] - stateOffsetMine; //number of states in current process
	size_t nStatesTot = stateOffset.back();
	
	//Separate energies by band:
	struct State
	{	double e, f, mfPrime; //energy, fermi filling, and -df/de
		vector3<> v; //band velocity
	};
	std::vector<State> state(nStatesMine);
	
	//N-fold division of k for error estimation:
	int nBlocks = 3;
	std::vector<int> omitIndex(h.nk); //which block to omit each k in
	std::vector<size_t> nkSel(nBlocks, h.nk); //k count per block after omissions
	Random::seed(0);
	for(int& omit: omitIndex)
	{	omit = Random::uniformInt(nBlocks);
		nkSel[omit]--;
	}
	mpiWorld->bcastData(omitIndex);
	mpiWorld->bcastData(nkSel);
	
	//Temperature loop:
	for(size_t iT=0; iT<Tcount; iT++)
	{	const double T = Tmin + iT*(Tmax-Tmin)/std::max(1, int(Tcount-1));
		const double invT = 1./T;
		
		std::vector<matrix3<>> rhoMat(nBlocks);
		std::vector<double> rho(nBlocks);
		std::vector<double> rhoRTA(nBlocks);
		std::vector<double> tauDrude(nBlocks);
		std::vector<double> tau(nBlocks);
		std::vector<double> vF(nBlocks);
		std::vector<double> gEf(nBlocks);
		
		for(int iBlock=0; iBlock<nBlocks; iBlock++)
		{	double omitWeight = h.nk*1./nkSel[iBlock];
			double vFsqSum = 0., weightSum = 0.;
			
			//Initialize states, velocity and -df/dE matrices:
			matrix V = zeroes(nStatesTot, 3);
			diagMatrix mfPrime(nStatesTot, 0.);
			State* s = state.data();
			size_t iState = stateOffsetMine;
			size_t nNZbound = 0; //upper bound number of non-zero e-ph matrix elements
			for(size_t ik=ikStart; ik<ikStop; ik++)
			{	const LindbladFile::Kpoint& kp = kpoint[ik-ikStart];
				for(int b=0; b<kp.nInner; b++)
				{	s->e = kp.E[kp.innerStart+b];
					double EminusMuByT = (s->e - dmu)*invT, fbar;
					fermi(EminusMuByT, s->f, fbar);
					mfPrime[iState] = s->mfPrime = (omitIndex[ik]==iBlock) ? 0. : invT * s->f * fbar; //project out k not in current block
					for(int iDir=0; iDir<3; iDir++)
						V.set(iState, iDir, s->v[iDir] = kp.P[iDir](b, kp.innerStart+b).real());
					//Fermi surface sums:
					vFsqSum += s->mfPrime * s->v.length_squared();
					weightSum += s->mfPrime;
					//e-ph matrix element count:
					size_t jk = -1;
					for(const LindbladFile::GePhEntry& g: kp.GePh)
					{	if(jk != g.jk) nNZbound += nInnerAll[g.jk];
						jk = g.jk;
					}
					s++;
					iState++;
				}
			}
			mpiWorld->allReduce(vFsqSum, MPIUtil::ReduceSum);
			mpiWorld->allReduce(weightSum, MPIUtil::ReduceSum);
			vF[iBlock] = sqrt(vFsqSum / weightSum);
			gEf[iBlock] = h.spinWeight*weightSum*omitWeight/h.nkTot;
			mpiWorld->allReduceData(V, MPIUtil::ReduceSum);
			mpiWorld->allReduceData(mfPrime, MPIUtil::ReduceSum);
			matrix mfPrimeV = mfPrime * V;
			
			//Construct BTE sparse matrix:
			logPrintf("\nConstructing BTE matrix: "); logFlush();
			TripletMatrix S(nStatesTot, nStatesTot, nNZbound, true);
			double prefacS = 2*M_PI*omitWeight/h.nkTot;
			double tauInvNum = 0.;
			for(size_t ik=ikStart; ik<ikStop; ik++)
			{	if(omitIndex[ik]!=iBlock)
				{	const LindbladFile::Kpoint& kp = kpoint[ik-ikStart];
					const State* s1 = state.data() + (stateOffset[ik]-stateOffsetMine);
					auto gStart = kp.GePh.begin();
					while(gStart != kp.GePh.end())
					{	size_t jk = gStart->jk;
						//Find range of g with same jk:
						auto gStop = gStart; gStop++;
						while((gStop != kp.GePh.end()) and (gStop->jk == jk)) gStop++;
						//Collect contributions:
						if(omitIndex[jk]!=iBlock)
						{	matrix Scontrib = zeroes(nInnerAll[ik], nInnerAll[jk]); //contributons to matrix Stilde in derivation
							complex* Sdata = Scontrib.data();
							for(auto g=gStart; g!=gStop; g++)
							{	double nPh = bose(invT*g->omegaPh); //phonon occupation
								for(const SparseEntry& e: g->G)
								{	double term = prefacS * e.val.norm(); //prefactors * |e-ph matrix element|^2 * energy conservation factor
									double nfbar_i = nPh + 1 - s1[e.i].f;
									double nf_j = nPh*(nPh+1)/nfbar_i; //nPh + fj by detailed balance
									Sdata[Scontrib.index(e.i,e.j)].real() -= term * (nfbar_i * s1[e.i].mfPrime);
									//Scattering time average:
									double fj = nf_j - nPh, mfjPrime = invT*fj*(1.-fj);
									tauInvNum += term * (nfbar_i * s1[e.i].mfPrime + nf_j * mfjPrime);
								}
							}
							for(int bj=0; bj<nInnerAll[jk]; bj++)
								for(int bi=0; bi<nInnerAll[ik]; bi++)
								{	if(Sdata->real())
									{	size_t iState = stateOffset[ik] + bi;
										size_t jState = stateOffset[jk] + bj;
										S.entries.push_back(TripletMatrix::Entry(iState,jState, Sdata->real()));
									}
									Sdata++;
								}
						}
						//Move to next jk set:
						gStart = gStop;
					}
				}
				//Print progress:
				if((ik-ikStart+1)%ikInterval==0) { logPrintf("%d%% ", int(round((ik-ikStart+1)*100./nkMine))); logFlush(); }
			}
			logPrintf("done.\n"); logFlush();
			mpiWorld->allReduce(tauInvNum, MPIUtil::ReduceSum);
			
			//Convert 'S' to matrix 'B' in the derivation:
			diagMatrix DiagSsum = S.DiagSumRows();
			for(size_t iState=stateOffsetMine; iState<stateOffsetMine+nStatesMine; iState++)
				S.entries.push_back(TripletMatrix::Entry(iState, iState, -0.5*DiagSsum[iState])); //factor of 0.5 due to implicit symmetry
			TripletMatrix& B = S; //call it B for clarity in the following
			
			//Relaxation time estimate:
			double Tt = (h.spinWeight*omitWeight/(3.*h.nkTot)) * trace(dagger(V) * mfPrimeV).real();
			double Gamma = (h.spinWeight*omitWeight/(3.*h.nkTot)) * trace(dagger(V) * (B * V)).real();
			rhoRTA[iBlock] = Omega*Gamma/(Tt*Tt);
			tauDrude[iBlock] = Tt / Gamma;
			
			//Calculate inv(B) * fPrimeV iteratively:
			double KmaxInv = 1e-6; //preconditioner regularizer
			for(size_t i=0; i<nStatesTot; i++)
			{	B.included[i] = mfPrime[i]; //only include non-zero points in solve subspace
				B.K[i] =  1./hypot(mfPrime[i], KmaxInv);
			}
			matrix invB_mfPrimeV = tauDrude[iBlock] * V; //initial guess
			B.applyInverse(mfPrimeV, invB_mfPrimeV);
			
			//Calculate conductivity tensor using Boltzmann equation:
			matrix3<> sigmaMat = (h.spinWeight*omitWeight/(h.nkTot*Omega)) * Matrix3(dagger(mfPrimeV) * invB_mfPrimeV);
			rhoMat[iBlock] = inv(sigmaMat);
			rho[iBlock] = (1./3) * trace(rhoMat[iBlock]);
			tau[iBlock] = weightSum / tauInvNum;
		}
		
		//Report:
		double rhoUnit = 1e-9*Ohm*meter;
		logPrintf("\nResults for T = %lg K:\n", T/Kelvin);
		reportResult(rhoMat, "Resistivity", rhoUnit, "nOhm-m");
		reportResult(rho, "Resistivity", rhoUnit, "nOhm-m");
		reportResult(rhoRTA, "ResistivityRTA", rhoUnit, "nOhm-m");
		reportResult(tauDrude, "tauDrude", fs, "fs");
		reportResult(tau, "tau", fs, "fs");
		reportResult(vF, "vF", 1, "");
		reportResult(gEf, "g(Ef)", 1, "");
		logPrintf("\n");
	}
	FeynWann::finalize();
}
