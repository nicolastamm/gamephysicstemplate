#include "DiffusionSimulator.h"
using namespace std;

Grid::Grid() :
	m(0),
	n(0)
{
	grid = std::vector<std::vector<Real>>();
}

void Grid::init(int m, int n)
{
	this->m = m;
	this->n = n;
	this->grid = std::vector<std::vector<Real>>(m, vector<Real>(n));
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			grid[i][j] = 0;
		}
	}
}

void Grid::fillRand()
{
	for (int i = 1; i < m - 1; i++)
	{
		for (int j = 1; j < n - 1; j++)
		{
			// fill with number [0, 0.999]
			grid[i][j] = rand() % 1000 / 1000.0;
		}
	}
}

int Grid::getM()
{
	return m;
}

int Grid::getN()
{
	return n;
}

std::vector<std::vector<Real>> Grid::getGrid()
{
	return grid;
}

Real Grid::read(int x, int y)
{
	return grid[x][y];
}

void Grid::write(int x, int y, Real val)
{
	grid[x][y] = val;
}


DiffusionSimulator::DiffusionSimulator()
{
	m_iTestCase = 0;
	m_vfMovableObjectPos = Vec3();
	m_vfMovableObjectFinalPos = Vec3();
	m_vfRotate = Vec3();
	// to be implemented
	T = new Grid();
	T->init(16, 16);
	T->fillRand();
	alpha = 1;
}

const char * DiffusionSimulator::getTestCasesStr() {
	return "Explicit_solver, Implicit_solver";
}

void DiffusionSimulator::reset() {
	m_mouse.x = m_mouse.y = 0;
	m_trackmouse.x = m_trackmouse.y = 0;
	m_oldtrackmouse.x = m_oldtrackmouse.y = 0;

}

void DiffusionSimulator::initUI(DrawingUtilitiesClass * DUC)
{
	this->DUC = DUC;
	// to be implemented
}

void DiffusionSimulator::notifyCaseChanged(int testCase)
{
	m_iTestCase = testCase;
	m_vfMovableObjectPos = Vec3(0, 0, 0);
	m_vfRotate = Vec3(0, 0, 0);
	//
	//to be implemented
	//
	switch (m_iTestCase)
	{
	case 0:
		cout << "Explicit solver!\n";
		break;
	case 1:
		cout << "Implicit solver!\n";
		break;
	default:
		cout << "Empty Test!\n";
		break;
	}
}

Grid* DiffusionSimulator::diffuseTemperatureExplicit(float timeStep) {//add your own parameters
	Grid* newT = new Grid();
	// to be implemented
	//make sure that the temperature in boundary cells stays zero
	int m = T->getM();
	int n = T->getN();
	newT->init(m, n);
	for (int i = 1; i < m - 1; i++) { // to avoid boundary values
		for (int j = 1; j < n - 1; j++) {
			Real t = T->read(i, j);
			Real t_xl = T->read(i, j - 1); // left
			Real t_xr = T->read(i, j + 1); // right
			Real t_yl = T->read(i - 1, j); // up
			Real t_yr = T->read(i + 1, j); // down
			// assume distance to be one here (need further check)
			Real x_partial = (t_xr - 2 * t + t_xl) / (1 * 1);
			Real y_partial = (t_yr - 2 * t + t_yl) / (1 * 1);
			Real t_new = t + timeStep * alpha * (x_partial + y_partial);
			newT->write(i, j, t_new);
		}
	}
	return newT;
}

void DiffusionSimulator::setupB(std::vector<Real>& b) {//add your own parameters
	// to be implemented
	//set vector B[sizeX*sizeY]
	int m = T->getM();
	int n = T->getN();

	for (int row = 0; row < m; row++)
	{
		for (int col = 0; col < n; col++)
			b.at(row * n + col) = T->read(row, col);	
	}

}

void fillT(Grid* T, vector<Real> x) {//add your own parameters
	// to be implemented
	//fill T with solved vector x
	//make sure that the temperature in boundary cells stays zero
	int m = T->getM();
	int n = T->getN();

	for (int row = 0; row < m; row++)
	{
		for (int col = 0; col < n; col++)
		{
			//First row is border, last row is border
			//First col is border, last col is border
			if (row == 0 || row == m - 1 || col == 0 || col == n - 1)
			{
				T->write(row, col, 0.0);
				continue;
			}
			T->write(row, col, x.at(row*n + col));
		}	
	}

}

void DiffusionSimulator::setupA(SparseMatrix<Real>& A, double factor) {//add your own parameters
	// to be implemented
	//setup Matrix A[sizeX*sizeY*sizeZ, sizeX*sizeY*sizeZ]
	// set with:  A.set_element( index1, index2 , value );
	// if needed, read with: A(index1, index2);
	// avoid zero rows in A -> set the diagonal value for boundary cells to 1.0
	int m = T->getM();
	int n = T->getN();

	for(int row = 0; row < m*n; row++)
	{
		//First m and last m rows are border.
		//Also, row % n == 0 and row % n == n-1 are borders too.
		//Set the diag value to one for all of these.
		if (row < m || row > pow(row, 2) - m
			|| row % n == 0 || row % n == n - 1)
		{
			A.set_element(row, row, 1.0);
			continue;
		}
		//Else --> Ai,i-1 = -F
		//         Ai,i   = 1 + 2F
		//         Ai,i+1 = -F
		A.set_element(row, row - 1, -factor);
		A.set_element(row, row, 1 + 2 * factor);
		A.set_element(row, row + 1, -factor);
	}
}


void DiffusionSimulator::diffuseTemperatureImplicit(float timeStep) {//add your own parameters
	// solve A T = b
	// to be implemented
	const int N = 16 * 16;//N = sizeX*sizeY*sizeZ
	SparseMatrix<Real> *A = new SparseMatrix<Real>(N);
	std::vector<Real> *b = new std::vector<Real>(N);

	setupA(*A, timeStep*alpha);
	setupB(*b);


	// perform solve
	Real pcg_target_residual = 1e-05;
	Real pcg_max_iterations = 1000;
	Real ret_pcg_residual = 1e10;
	int  ret_pcg_iterations = -1;

	SparsePCGSolver<Real> solver;
	solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97, 0.25);

	std::vector<Real> x(N);
	for (int j = 0; j < N; ++j) { x[j] = 0.; }

	// preconditioners: 0 off, 1 diagonal, 2 incomplete cholesky
	solver.solve(*A, *b, x, ret_pcg_residual, ret_pcg_iterations, 2);
	// x contains the new temperature values
	fillT(T,x);//copy x to T
}

void DiffusionSimulator::onGridSizeChange(int m, int n)
{
	// TODO: implement this
}



void DiffusionSimulator::simulateTimestep(float timeStep)
{
	// to be implemented
	// update current setup for each frame
	switch (m_iTestCase)
	{
	case 0:
		T = diffuseTemperatureExplicit(timeStep);
		break;
	case 1:
		diffuseTemperatureImplicit(timeStep);
		break;
	}
}

void DiffusionSimulator::drawObjects()
{
	// to be implemented
	//visualization
	int m = T->getM();
	int n = T->getN();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			Real t = T->read(i, j);
			DUC->setUpLighting(Vec3(0, 0, 0), 0.4 * Vec3(1, 1, 1), 100., Vec3(t, 0, 1 - t));
			DUC->drawSphere(Vec3(i, j, 0), Vec3(1, 1, 1));
		}
	}
}


void DiffusionSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext)
{
	drawObjects();
}

void DiffusionSimulator::onClick(int x, int y)
{
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

void DiffusionSimulator::onMouse(int x, int y)
{
	m_oldtrackmouse.x = x;
	m_oldtrackmouse.y = y;
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}