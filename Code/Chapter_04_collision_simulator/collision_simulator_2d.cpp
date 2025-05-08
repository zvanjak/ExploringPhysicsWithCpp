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
	struct Ball2D
	{
	private:
		double _mass;
		double _radius;
		std::string _color;
		Point2Cartesian _position;
		Vector2Cartesian _velocity;

	public:
		Ball2D(double mass, double radius, std::string color, const Point2Cartesian& position,
			const Vector2Cartesian& velocity)
			: _mass(mass), _radius(radius), _color(color), _position(position), _velocity(velocity)
		{
		}

		double  Mass() const { return _mass; }
		double& Mass() { return _mass; }

		double  Rad() const { return _radius; }
		double& Rad() { return _radius; }

		std::string Color() const { return _color; }
		std::string& Color() { return _color; }

		Point2Cartesian  Pos() const { return _position; }
		Point2Cartesian& Pos() { return _position; }

		Vector2Cartesian  V() const { return _velocity; }
		Vector2Cartesian& V() { return _velocity; }
	};

	struct Container2D
	{
		double _width;
		double _height;

		std::vector<Ball2D> _balls;

		Container2D() : _width(1000), _height(1000) {}
		Container2D(double width, double height) : _width(width), _height(height) {}

		// add body
		void AddBall(const Ball2D& body)
		{
			_balls.push_back(body);
		}

		// get body
		Ball2D& Ball(int i)
		{
			return _balls[i];
		}

		// check if ball is out of bounds and handle it
		void CheckAndHandleOutOfBounds(int ballIndex)
		{
			Ball2D& ball = _balls[ballIndex];

			// left wall collision
			if (ball.Pos().X() < ball.Rad() && ball.V().X() < 0)
			{
				ball.Pos().X() = ball.Rad() + (ball.Rad() - ball.Pos().X()); // Get back to box!
				ball.V().X() *= -1;
			}
			// right wall collision
			if (ball.Pos().X() > _width - ball.Rad() && ball.V().X() > 0)
			{
				ball.Pos().X() -= (ball.Pos().X() + ball.Rad()) - _width;
				ball.V().X() *= -1;
			}
			// bottom wall collision
			if (ball.Pos().Y() < ball.Rad() && ball.V().Y() < 0)
			{
				ball.Pos().Y() = ball.Rad() + (ball.Rad() - ball.Pos().Y());
				ball.V().Y() *= -1;
			}
			// top wall collision
			if (ball.Pos().Y() > _height - ball.Rad() && ball.V().Y() > 0)
			{
				ball.Pos().Y() -= (ball.Pos().Y() + ball.Rad()) - _height;
				ball.V().Y() *= -1;
			}
		}
	};

	class CollisionSimulator2D
	{
		Container2D _box;

	public:
		CollisionSimulator2D() {}
		CollisionSimulator2D(const Container2D& box) : _box(box) {}

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
						HandleCollision(m, n, dt);
				}
			}
		}

		void HandleCollision(int m, int n, double dt)
		{
			Ball2D& ball1 = _box._balls[m];
			Ball2D& ball2 = _box._balls[n];

			// calculating point where they were before collision
			// (if there was collision with the box wall, then calc.pos. will be outside the box
			// but it doesn't matter, since we need only direction, ie. velocity, to calculate exact collision point)
			Pnt2Cart x10 = ball1.Pos() - ball1.V() * dt;
			Pnt2Cart x20 = ball2.Pos() - ball2.V() * dt;

			Vec2Cart dx0(x10, x20);
			Vec2Cart dv(ball2.V() - ball1.V());

			// first, we have to calculate exact moment of collision, and balls positions then
			double A = dv * dv;
			double B = 2 * dx0 * dv;
			double C = dx0 * dx0 - POW2(ball2.Rad() + ball1.Rad());

			double t1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
			double t2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);

			double tCollision = t1 < t2 ? t1 : t2;
			//double tReverseBalls = tCollision - dt;    // with respect to CURRENT balls position

			// calculating position of balls at the point of collision (moving them backwards)
			Pnt2Cart x1 = ball1.Pos() + (tCollision - dt) * ball1.V();
			Pnt2Cart x2 = ball2.Pos() + (tCollision - dt) * ball2.V();

			// https://en.wikipedia.org/wiki/Elastic_collision - calculating new velocities after collision
			double   m1 = ball1.Mass(), m2 = ball2.Mass();

			Vec2Cart v1 = ball1.V(), v2 = ball2.V();

			Vec2Cart v1_v2 = v1 - v2;
			Vec2Cart x1_x2(x2, x1);

			Vec2Cart v1_new = v1 - 2 * m2 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec2Cart(x2, x1);
			Vec2Cart v2_new = v2 - 2 * m1 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec2Cart(x1, x2);

			ball1.V() = v1_new;
			ball2.V() = v2_new;

			// adjusting new ball positions
			ball1.Pos() = x1 + ball1.V() * (dt - tCollision);
			ball2.Pos() = x2 + ball2.V() * (dt - tCollision);
		}

		std::vector<std::vector<Pnt2Cart>> Simulate(int numSteps, double timeStep)
		{
			int numBalls = _box._balls.size();

			std::vector<std::vector<Pnt2Cart>> ballPositions(numBalls);

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

		void Serialize(std::string fileName, std::vector<std::vector<Pnt2Cart>> ballPositions)
		{
			std::ofstream file(fileName);
			if (file.is_open())
			{
				file << "PARTICLE_SIMULATION_DATA_2D" << std::endl;
				file << "NumBalls: " << _box._balls.size() << std::endl;

				int count = 0;
				for (const auto& ball : _box._balls)
				{
					file << "Ball_" << ++count << " " << ball.Color() << " " << ball.Rad() << std::endl;
				}

				file << "NumSteps: " << ballPositions[0].size() << std::endl;

				for (int i = 0; i < ballPositions[0].size(); i++)
				{
					file << "Step " << i << " 0.1" << std::endl;
					for (int j = 0; j < _box._balls.size(); j++)
					{
						const Ball2D& ball = _box._balls[j];
						file << j << " " << ballPositions[j][i].X() << " " << ballPositions[j][i].Y() << "\n";
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

double RandomUniform(double min, double max)
{
	return min + (max - min) * ((double)rand() / RAND_MAX);
}

void Example4_collision_calculator_2D()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                 EXAMPLE 4 - collision calculator              ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// example 1
	Container2D box1(1000, 800);

	box1.AddBall(Ball2D(1.0, 10.0, "red", Pnt2Cart(100, 150), Vec2Cart(3, 1.5)));
	box1.AddBall(Ball2D(2.0, 20.0, "blue", Pnt2Cart(200, 500), Vec2Cart(-2, -3)));
	box1.AddBall(Ball2D(1.0, 10.0, "green", Pnt2Cart(300, 150), Vec2Cart(1, -2)));
	box1.AddBall(Ball2D(3.0, 20.0, "yellow", Pnt2Cart(100, 400), Vec2Cart(-1, 1)));
	box1.AddBall(Ball2D(1.0, 10.0, "purple", Pnt2Cart(200, 600), Vec2Cart(1, -0.6)));

	// example 2 - 100 random balls
	Container2D box2(1000, 800);
	for (int i = 0; i < 50; i++)
	{
		double mass = RandomUniform(1.0, 10.0);
		double radius = RandomUniform(2.0, 10.0);
		std::string color = "Red";
		Point2Cartesian position(RandomUniform(radius, box2._width - radius), RandomUniform(radius, box2._height - radius));
		Vector2Cartesian velocity(RandomUniform(-5.0, 5.0), RandomUniform(-5.0, 5.0));

		box2.AddBall(Ball2D(mass, radius, color, position, velocity));
	}
	for (int i = 0; i < 50; i++)
	{
		double mass = RandomUniform(1.0, 10.0);
		double radius = RandomUniform(2.0, 10.0);
		std::string color = "Blue";
		Point2Cartesian position(RandomUniform(radius, box2._width - radius), RandomUniform(radius, box2._height - radius));
		Vector2Cartesian velocity(RandomUniform(-5.0, 5.0), RandomUniform(-5.0, 5.0));

		box2.AddBall(Ball2D(mass, radius, color, position, velocity));
	}

	// create simulator
	CollisionSimulator2D simulator(box2);

	// simulate
	int numSteps = 1000;
	auto ballPositions = simulator.Simulate(numSteps, 1.0);

	simulator.Serialize(MML_PATH_ResultFiles + "SimData2.txt", ballPositions);

	Visualizer::VisualizeParticleSimulation2D("SimData2.txt");

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