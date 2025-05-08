#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorTypes.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"
#endif

using namespace MML;

// Calculating deviation and changes for ever decreasing dT
//   naci onu tocku nakon koje smanjivanje dT ne vodi drugacijim rezultatima

namespace CollisionSimulator
{
	struct Ball3D
	{
	private:
		double _mass;
		double _radius;
		std::string _color;
		Pnt3Cart _position;
		Vec3Cart _velocity;

	public:
		Ball3D(double mass, double radius, std::string color, const Pnt3Cart& position,
			const Vec3Cart& velocity)
			: _mass(mass), _radius(radius), _color(color), _position(position), _velocity(velocity)
		{
		}

		double  Mass() const { return _mass; }
		double& Mass() { return _mass; }

		double  Rad() const { return _radius; }
		double& Rad() { return _radius; }

		std::string Color() const { return _color; }
		std::string& Color() { return _color; }

		Pnt3Cart  Pos() const { return _position; }
		Pnt3Cart& Pos() { return _position; }

		Vec3Cart  V() const { return _velocity; }
		Vec3Cart& V() { return _velocity; }
	};

	struct Container3D
	{
		double _width;
		double _height;
		double _depth;

		std::vector<Ball3D> _balls;

		Container3D() : _width(1000), _height(1000), _depth(1000) {}
		Container3D(double width, double height, double depth) : _width(width), _height(height), _depth(depth) {}

		// add body
		void AddBall(const Ball3D& body)
		{
			_balls.push_back(body);
		}

		// get body
		Ball3D& Ball(int i)
		{
			return _balls[i];
		}

		// check if ball is out of bounds and handle it

		void CheckAndHandleOutOfBounds(int ballIndex)
		{
			Ball3D& ball = _balls[ballIndex];

			if (ball.Pos().X() < ball.Rad() && ball.V().X() < 0)		// checking if ball is out of bounds
			{
				ball.Pos().X() = ball.Rad() + (ball.Rad() - ball.Pos().X());		// moving it back to box
				ball.V().X() *= -1;
			}

			if (ball.Pos().X() > _width - ball.Rad() && ball.V().X() > 0)
			{
				ball.Pos().X() -= (ball.Pos().X() + ball.Rad()) - _width;
				ball.V().X() *= -1;
			}

			if (ball.Pos().Y() < ball.Rad() && ball.V().Y() < 0)
			{
				ball.Pos().Y() = ball.Rad() + (ball.Rad() - ball.Pos().Y());
				ball.V().Y() *= -1;
			}

			if (ball.Pos().Y() > _height - ball.Rad() && ball.V().Y() > 0)
			{
				ball.Pos().Y() -= (ball.Pos().Y() + ball.Rad()) - _height;
				ball.V().Y() *= -1;
			}

			if (ball.Pos().Z() < ball.Rad() && ball.V().Z() < 0)
			{
				ball.Pos().Z() = ball.Rad() + (ball.Rad() - ball.Pos().Z());
				ball.V().Z() *= -1;
			}

			if (ball.Pos().Z() > _depth - ball.Rad() && ball.V().Z() > 0)
			{
				ball.Pos().Z() -= (ball.Pos().Z() + ball.Rad()) - _depth;
				ball.V().Z() *= -1;
			}
		}
	};

	class CollisionSimulator3D
	{
		Container3D _box;

	public:
		CollisionSimulator3D() {}
		CollisionSimulator3D(const Container3D& box) : _box(box) {}

		double DistBalls(int i, int j)
		{
			return _box._balls[i].Pos().Dist(_box._balls[j].Pos());
		}

		bool HasBallsCollided(int i, int j)
		{
			// ako je udaljenost izmedju njihovih centara manja od zbroja radijusa
			if (DistBalls(i, j) < _box._balls[i].Rad() + _box._balls[j].Rad())
				return true;
			else
				return false;
		}

		std::vector<std::vector<Pnt3Cart>> Simulate(int numSteps, double timeStep)
		{
			int numBalls = _box._balls.size();

			std::vector<std::vector<Pnt3Cart>> ballPositions(numBalls);

			for (int i = 0; i < numSteps; i++)
			{
				double dt = timeStep; // time step

				// save positions
				for (int j = 0; j < numBalls; j++)
					ballPositions[j].push_back(_box._balls[j].Pos());

				// simulate one step
				SimulateOneStep(dt);
			}

			return ballPositions;
		}

		void SimulateOneStep(double dt)
		{
			// check collisions at start configuration!!!
			// that is also an end configuration for previous step
			int NumBalls = _box._balls.size();

			// first, update all ball's positions, and handle out of bounds
			for (int i = 0; i < NumBalls; i++)
			{
				_box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;

				_box.CheckAndHandleOutOfBounds(i);
			}

			// check for collisions, and handle it if there is one
			for (int m = 0; m < NumBalls - 1; m++) {
				for (int n = m + 1; n < NumBalls; n++)
				{
					if (HasBallsCollided(m, n))
					{
						std::cout << "Collision detected between balls " << m << " and " << n << std::endl;

						Ball3D& ball1 = _box._balls[m];
						Ball3D& ball2 = _box._balls[n];

						// calculating point where they were before collision
						// (if there was collision with the box wall, then calc.pos. will be outside the box
						// but it doesn't matter, since we need only direction, ie. velocity, to calculate exact collision point)
						Pnt3Cart x10 = ball1.Pos() - ball1.V() * dt;
						Pnt3Cart x20 = ball2.Pos() - ball2.V() * dt;

						Vec3Cart dx0(x10, x20);
						Vec3Cart dv(ball2.V() - ball1.V());

						double A = dv * dv;
						double B = 2 * dx0 * dv;
						double C = dx0 * dx0 - POW2(ball2.Rad() + ball1.Rad());

						double t1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
						double t2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);

						double tCollision = t1 < t2 ? t1 : t2;
						//double tReverseBalls = tCollision - dt;    // with respect to CURRENT balls position

						// calculating position of balls at the point of collision (moving them backwards)
						Pnt3Cart x1 = ball1.Pos() + (tCollision - dt) * ball1.V();
						Pnt3Cart x2 = ball2.Pos() + (tCollision - dt) * ball2.V();

						// https://en.wikipedia.org/wiki/Elastic_collision - calculating new velocities after collision
						double   m1 = ball1.Mass(), m2 = ball2.Mass();

						Vec3Cart v1 = ball1.V(), v2 = ball2.V();

						Vec3Cart v1_v2 = v1 - v2;
						Vec3Cart x1_x2(x2, x1);

						Vec3Cart v1_new = v1 - 2 * m2 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec3Cart(x2, x1);
						Vec3Cart v2_new = v2 - 2 * m1 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec3Cart(x1, x2);

						ball1.V() = v1_new;
						ball2.V() = v2_new;

						// adjusting new ball positions
						ball1.Pos() = x1 + ball1.V() * (dt - tCollision);
						ball2.Pos() = x2 + ball2.V() * (dt - tCollision);
					}
				}
			}
		}

		void Serialize(std::string fileName, std::vector<std::vector<Pnt3Cart>> ballPositions)
		{
			std::ofstream file(fileName);
			if (file.is_open())
			{
				file << "PARTICLE_SIMULATION_DATA_3D" << std::endl;
				file << "NumBalls: " << _box._balls.size() << std::endl;

				for (const auto& ball : _box._balls)
				{
					file << "Ball_" << " " << ball.Color() << " " << ball.Rad() << std::endl;
				}

				file << "NumSteps: " << ballPositions[0].size() << std::endl;

				for (int i = 0; i < ballPositions[0].size(); i++)
				{
					file << "Step " << i << " 0.1" << std::endl;
					for (int j = 0; j < _box._balls.size(); j++)
					{
						const Ball3D& ball = _box._balls[j];
						file << j << " " << ballPositions[j][i].X() << " " << ballPositions[j][i].Y() << " " << ballPositions[j][i].Z() << "\n";
					}
				}
				file.close();
			}
			else
			{
				std::cerr << "Unable to open file";
			}
		}
	};
}

using namespace CollisionSimulator;


void Example4_collision_calculator_3D()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                 EXAMPLE 4 - collision calculator 3D           ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// example 1
	Container3D box1(1000, 1000, 1000);

	box1.AddBall(Ball3D(1.0, 10.0, "red", Pnt3Cart(100, 150, 50), Vec3Cart(3, 1.5, 4)));
	box1.AddBall(Ball3D(2.0, 20.0, "blue", Pnt3Cart(200, 500, 200), Vec3Cart(-2, -3, -6)));
	box1.AddBall(Ball3D(1.0, 10.0, "green", Pnt3Cart(300, 150, 300), Vec3Cart(1, -2, 3)));
	box1.AddBall(Ball3D(3.0, 20.0, "yellow", Pnt3Cart(100, 400, 200), Vec3Cart(-1, 1, 3)));
	box1.AddBall(Ball3D(1.0, 10.0, "purple", Pnt3Cart(200, 600, 500), Vec3Cart(1, -0.6, -2)));

	// example 2 - 100 random balls
	Container3D box2(1000, 1000, 1000);

	std::vector<std::string> vecColors(100);
	std::vector<double> vecRad(100);

	for (int i = 0; i < 50; i++)
	{
		double mass = Random::UniformReal(1.0, 10.0);
		double radius = Random::UniformReal(2.0, 10.0);
		vecRad[i] = radius;
		vecColors[i] = "Red";
		std::string color = "Red";
		Pnt3Cart position(Random::UniformReal(radius, box2._width - radius), Random::UniformReal(radius, box2._height - radius), Random::UniformReal(radius, box2._height - radius));
		Vec3Cart velocity(Random::UniformReal(-5.0, 5.0), Random::UniformReal(-5.0, 5.0), Random::UniformReal(-5.0, 5.0));

		box2.AddBall(Ball3D(mass, radius, color, position, velocity));
	}
	for (int i = 0; i < 50; i++)
	{
		double mass = Random::UniformReal(1.0, 10.0);
		double radius = Random::UniformReal(2.0, 10.0);
		vecRad[50 + i] = radius;
		std::string color = "Blue";
		vecColors[50 + i] = color;
		Pnt3Cart position(Random::UniformReal(radius, box2._width - radius), Random::UniformReal(radius, box2._height - radius), Random::UniformReal(radius, box2._height - radius));
		Vec3Cart velocity(Random::UniformReal(-5.0, 5.0), Random::UniformReal(-5.0, 5.0), Random::UniformReal(-5.0, 5.0));

		box2.AddBall(Ball3D(mass, radius, color, position, velocity));
	}

	// create simulator
	CollisionSimulator3D simulator(box2);

	// simulate
	int numSteps = 100;
	auto ballPositions = simulator.Simulate(numSteps, 1.0);

	Serializer::SaveParticleSimulation3D(MML_PATH_ResultFiles + "SimData3D.txt", box2._balls.size(), ballPositions,
		vecColors, vecRad);

	//simulator.Serialize(MML_PATH_ResultFiles + "SimData3.txt", ballPositions);

	Visualizer::VisualizeParticleSimulation3D("SimData3D.txt");

	// print results

	//for (int i = 0; i < numSteps; i++)
	//{
	//	std::cout << "Step " << i << ": " << std::endl;
	//	for (int j = 0; j < box._balls.size(); j++)
	//	{
	//		std::cout << "Ball " << j << ": " << ballPositions[j][i].X() << ", " << ballPositions[j][i].Y() << std::endl;
	//	}
	//}
}