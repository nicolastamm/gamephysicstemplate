#ifndef DIFFUSIONSIMULATOR_h
#define DIFFUSIONSIMULATOR_h

#include "Simulator.h"
#include "vectorbase.h"
#include "util/FFmpeg.h"
#include "pcgsolver.h"

//impement your own grid class for saving grid data
class Grid {
public:
	// Construtors
	Grid();
	void init(int m, int n);
	void fillRand();

	int getM();
	int getN();
	std::vector<std::vector<Real>> getGrid();
	Real read(int x, int y);
	void write(int x, int y, Real val);

private:
	// Attributes
	std::vector<std::vector<Real>> grid;
	int m;
	int n;
};



class DiffusionSimulator :public Simulator {
public:
	// Construtors
	DiffusionSimulator();

	// Functions
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
	void simulateTimestep(float timeStep);
	void externalForcesCalculations(float timeElapsed) {};
	void onClick(int x, int y);
	void onMouse(int x, int y);
	// Specific Functions
	void drawObjects();
	Grid* diffuseTemperatureExplicit(float timeStep);
	void diffuseTemperatureImplicit(float timeStep);
	void onGridSizeChange(int m, int n);

private:
	// Attributes
	Vec3  m_vfMovableObjectPos;
	Vec3  m_vfMovableObjectFinalPos;
	Vec3  m_vfRotate;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
	Grid *T; //save results of every time step
	Real alpha;
	// methods
	void setupB(std::vector<Real>& b);
	void setupA(SparseMatrix<Real>& A, double factor);
};

#endif