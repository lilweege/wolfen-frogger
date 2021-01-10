#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <assert.h>

/*
* https://lodev.org/cgtutor/raycasting.html#Introduction
* https://brilliant.org/wiki/dot-product-distance-between-point-and-a-line/
* https://thecodingtrain.com/CodingChallenges/146-rendering-ray-casting.html
* https://github.com/OneLoneCoder/videos/blob/master/OneLoneCoder_ComandLineFPS_2.cpp
*/

static const olc::vd2d NOVD2D = { DBL_MIN, DBL_MIN };

static bool wrap(double& x, double m)
{
	if (x > m)
	{
		x -= m;
		return true;
	}

	if (x < 0)
	{
		x += m;
		return true;
	}

	return false;
}

static olc::vd2d lineLine(const olc::vd2d& p1, const olc::vd2d& p2, const olc::vd2d& p3, const olc::vd2d& p4)
{
	const double
		x43 = p4.x - p3.x,
		x13 = p1.x - p3.x,
		x21 = p2.x - p1.x,
		y43 = p4.y - p3.y,
		y13 = p1.y - p3.y,
		y21 = p2.y - p1.y;

	const double d = y43 * x21 - x43 * y21;
	const double uA = (x43 * y13 - y43 * x13) / d;
	const double uB = (x21 * y13 - y21 * x13) / d;

	if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1)
		return { p1.x + uA * x21, p1.y + uA * y21 };

	return NOVD2D;
}

static bool pointPoly(const olc::vd2d& p, const olc::vd2d& q, const std::vector<std::pair<olc::vd2d, olc::vd2d>>& polygon)
{
	bool collision = false;
	for (auto [a, b] : polygon)
	{
		a += q, b += q;
		if (((a.y >= p.y && b.y < p.y) || (a.y < p.y && b.y >= p.y)) &&
			(p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x))
			collision = !collision;
	}
	
	return collision;
}

class olcWolfenFrogger : public olc::PixelGameEngine
{
public:
	olcWolfenFrogger()
	{
		sAppName = "Wolfen-Frogger";
	}

private:

	struct Entity2D
	{
		olc::vd2d position;
		olc::vd2d velocity;

		double movementSpeed;
		double rotationSpeed;

		// velocity does not necessarily imply direction
		double heading = 0;
		olc::vd2d direction = { 1, 0 };


		Entity2D(const olc::vd2d& pos, const olc::vd2d& vel = { 0, 0 }, double mov = 0, double rot = 0)
			: position(pos), velocity(vel), movementSpeed(mov), rotationSpeed(rot) {}

		void setPos(const olc::vd2d& pos)
		{
			position = pos;
		}

		void setHead(double head)
		{
			heading = head;
			direction.x = cos(heading);
			direction.y = sin(heading);
		}


		// formally move and rotate

		void walk(float step)
		{
			position += velocity * movementSpeed * step;
		}

		void turn(int dir, float step)
		{
			heading += dir * rotationSpeed * step;
			direction.x = cos(heading);
			direction.y = sin(heading);
		}
	};

	struct WorldObject : public Entity2D
	{
		int numFaces;
		std::vector<std::pair<olc::vd2d, olc::vd2d>> polygon;
		//std::vector<olc::Decal*> textures;
		olc::Pixel color;

		//std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> faceEdges;
		//std::vector<bool> faceVisible;
		
		// polygon points should be in CCW order for texture to be correct
		// texture will be flipped if face is viewed from the back side
		WorldObject(const olc::vd2d& pos,
			const std::vector<olc::vd2d>& points,
			const olc::Pixel& col = olc::WHITE/*,
			const std::vector<olc::Decal*>& decals = {}*/,
			const olc::vd2d& vel = { 0, 0 },
			double mov = 0, double rot = 0)
			: Entity2D(pos, vel, mov, rot), color(col), numFaces(points.size())
		{
			// don't make dumb objects pls
			assert(numFaces > 0);

			polygon.reserve(numFaces);
			//textures.reserve(numFaces);
			//int numDecals = decals.size();
			for (int i = 0; i < numFaces; ++i)
			{
				polygon.push_back(std::make_pair(points[i], points[i == numFaces - 1 ? 0 : i + 1]));
				//textures.push_back(i >= numDecals ? nullptr : decals[i]);
			}
			//faceEdges.resize(numFaces);
			//faceVisible.resize(numFaces);
		}
	};


	struct Car : public WorldObject
	{
		// MUST BE PASSED: pos, dim -> points, color
		// position, size, 
		// initial velocity, initial angle, movementSpeed, rotationSpeed
		Car(const olc::vd2d& pos, const olc::vd2d& dim,
			const olc::vd2d& vel, double speed,
			const olc::Pixel& col = olc::WHITE/*,
			const std::vector<olc::Decal*>& decals = {}*/)
			: WorldObject(pos,
				{	
					olc::vd2d(    0,     0),
					olc::vd2d(    0, dim.y),
					olc::vd2d(dim.x, dim.y),
					olc::vd2d(dim.x,     0)
				}, col/*, decals*/, vel, speed)
		{
		}
	};


	struct Player : public Entity2D
	{
		double fieldOfView;
		double rayLength;
		double planeDist;
		double planeLength;

		struct Ray
		{
			olc::vd2d target = NOVD2D;
			olc::vd2d hit = NOVD2D;
			double length = DBL_MAX;

			// color information stored per ray
			olc::Pixel col = olc::WHITE;
			// when textures are added some sort
			// of z-buffer may be necessary
		};
		int numRays;
		std::vector<Ray> rays;

		Player(const olc::vi2d& worldSize, double worldScale,
			const olc::vd2d& pos = { 0, 0 },
			double head = 0, double FOV = M_PI_2)
			: Entity2D(pos, { 0, 0 }, worldSize.x / 4.0 / worldScale, 2)
		{
			numRays = worldSize.x / 2; // arbitrary
			rays.resize(numRays);

			planeDist = worldSize.x / worldScale;
			fieldOfView = FOV;
			rayLength = planeDist / cos(fieldOfView / 2.0);
			planeLength = 2 * sqrt(rayLength * rayLength - planeDist * planeDist);

			position = pos;
			setHead(head);
			recalcRays();
		}

		void update(const olcWolfenFrogger* game, double timestep)
		{
			// don't be dum
			assert(game != nullptr);

			bool
				wPress = game->GetKey(olc::Key::W).bHeld,
				aPress = game->GetKey(olc::Key::A).bHeld,
				sPress = game->GetKey(olc::Key::S).bHeld,
				dPress = game->GetKey(olc::Key::D).bHeld,
				qPress = game->GetKey(olc::Key::Q).bHeld,
				ePress = game->GetKey(olc::Key::E).bHeld;

			int rotDir = ePress - qPress;
			if (rotDir)
				turn(rotDir, timestep);
			
			velocity = (
				+ direction * wPress
				- direction * sPress
				+ direction.perp() * dPress
				- direction.perp() * aPress
			);
			// avoid norm div by 0
			if (velocity.x && velocity.y)
			{
				velocity = velocity.norm();
				walk(timestep);
			}

			if (wPress || aPress || sPress || dPress || qPress || ePress)
				recalcRays();
			castRays(game->objs);
		}


		void recalcRays()
		{
			for (int i = 0; i < numRays; ++i)
				rays[i].target =
					position +
					direction * planeDist +
					direction.perp() * planeLength * (i / double(numRays) - 0.5);
		}

		void castRays(const std::vector<WorldObject*>& objs)
		{
			//for (WorldObject* obj : objs)
			//	std::fill(obj->faceVisible.begin(), obj->faceVisible.end(), false);

			const olc::vd2d& p3 = position;
			for (Ray& ray : rays)
			{
				const olc::vd2d& p4 = ray.target;
				ray.length = DBL_MAX;
				ray.hit = NOVD2D;
				ray.col = olc::WHITE;

				//WorldObject* hit = nullptr;
				//int face;

				for (WorldObject* obj : objs)
					for (int f = 0; f < obj->numFaces; ++f)
					{
						olc::vd2d p1 = obj->position + obj->polygon[f].first;
						olc::vd2d p2 = obj->position + obj->polygon[f].second;

						olc::vd2d intersect = lineLine(p1, p2, p3, p4);

						if (intersect != NOVD2D)
						{
							// projected perpendicular distance
							// euclidian results in fisheye
							double distance = (intersect - position).dot(direction);

							if (ray.length > distance)
							{
								ray.length = distance;
								ray.hit = intersect;
								ray.col = obj->color;

								//hit = obj;
								//face = f;
							}
						}
					}

				//if (hit != nullptr)
				//	hit->faceVisible[face] = true;
			}

			/*
			for (WorldObject* obj : objs)
				for (int f = 0; f < obj->numFaces; ++f)
				{
					if (!obj->faceVisible[f])
						continue;

					// left and right edges

					olc::vd2d pl = obj->position + obj->polygon[f].first - position;
					olc::vd2d pr = obj->position + obj->polygon[f].second - position;

					double dl = pl.dot(direction);
					double dr = pr.dot(direction);

					auto& [el, er] = obj->faceEdges[f];
					// projected distance -> height
					el.second = std::min(std::max(0.0, dl), planeDist);
					er.second = std::min(std::max(0.0, dr), planeDist);
					// "index position" -> x pos
					el.first = pl.dot(direction.perp()) / el.second * numRays / 2 + numRays / 2;
					er.first = pr.dot(direction.perp()) / er.second * numRays / 2 + numRays / 2;
				}
			*/
		}
	};

private:
	double worldScale = 10;
	double worldWidth = ScreenWidth() / worldScale;
	double worldHeight = ScreenHeight() / worldScale;
	double minimapScale = worldScale * 0.25;

	Player* frog;
	olc::vd2d startPos;
	double startAngle;

	Car* car1;
	Car* car2;
	Car* car3;
	std::vector<WorldObject*> objs;

	olc::Pixel* bgBuffer;
	void Clear()
	{
		int pixels = GetDrawTargetWidth() * GetDrawTargetHeight();
		olc::Pixel* m = GetDrawTarget()->GetData();
		for (int i = 0; i < pixels; ++i)
			m[i] = bgBuffer[i];
	}

public:
	bool OnUserCreate() override
	{
		//precompute sky and floor background
		bgBuffer = (olc::Pixel*) malloc(ScreenWidth() * ScreenHeight() * sizeof(olc::Pixel));

		// sky
		const olc::Pixel skyColor = olc::Pixel(127, 255, 255);
		for (int x = 0; x < ScreenWidth(); ++x)
			for (int y = 0; y < ScreenHeight() / 2; ++y)
				bgBuffer[x + y * ScreenWidth()] = skyColor;

		// floor
		const olc::Pixel floorColor = olc::GREY;
		for (int y = ScreenHeight() / 2; y < ScreenHeight(); ++y)
		{
			float bright = float(1 - sqrt(ScreenHeight() - y) / 10);
			for (int x = 0; x < ScreenWidth(); ++x)
				bgBuffer[x + y * ScreenWidth()] = floorColor * bright;
		}



		/*init game state*/

		// start bottom middle facing upwards
		startPos = { worldWidth * 0.5, worldHeight - 2 / worldScale };
		startAngle = -M_PI_2;
		frog = new Player({ ScreenWidth(), ScreenHeight() }, worldScale, startPos, startAngle);

		// extra wall border
		double
			x = worldWidth - 1 / worldScale,
			y = worldHeight - 1 / worldScale;
		objs.push_back(
			new WorldObject({ 0, 0 }, {
				{ 0, 0 },
				{ x, 0 },
				{ x, y },
				{ 0, y }
			})
		);

		olc::vd2d carDim = { worldWidth / 15.0, worldHeight / 15.0 };
		car1 = new Car({ 0, 3 }, carDim, { 1, 0 }, 5.0, olc::RED);
		objs.push_back(car1);

		car2 = new Car({ 0, 7 }, carDim, { -1, 0 }, 10.0, olc::GREEN);
		objs.push_back(car2);

		car3 = new Car({ 0, 11 }, carDim, { 1, 0 }, 3.0, olc::BLUE);
		objs.push_back(car3);


		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		/*logic stuff*/

		// update stuff
		for (WorldObject* obj : objs)
		{
			obj->walk(fElapsedTime);

			// fmod is suck
			wrap(obj->position.x, worldWidth);
			wrap(obj->position.y, worldHeight);
		}

		// player update
		frog->update(this, fElapsedTime);

		// collide
		if (pointPoly(frog->position, car1->position, car1->polygon) ||
			pointPoly(frog->position, car3->position, car3->polygon) ||
			pointPoly(frog->position, car2->position, car2->polygon))
		{
			// reset
			frog->setPos(startPos);
			frog->setHead(startAngle);
			frog->recalcRays();
		}


		/*draw stuff*/
		Clear();

		// cast rays
		const int mY = ScreenHeight() / 2;
		for (int x = 0; x < ScreenWidth(); ++x)
		{
			int i = int(x * frog->numRays / float(ScreenWidth()));
			const Player::Ray& ray = frog->rays[i];
			if (ray.hit == NOVD2D)
				continue;

			const double dist = std::min(std::max(1.0, ray.length), frog->planeDist);
			const int h = int(ScreenHeight() / dist);

			float c = float(ScreenHeight() / frog->planeDist);
			float brightness = sqrt((h - c) / (float(ScreenHeight()) - c));

			DrawLine(x, mY - h / 2, x, mY + h / 2, ray.col * brightness);
 		}


		// textures ???
		/*
		// FACES DRAWN WITH ARBITRARY Z-ORDER
		for (const WorldObject* obj : objs)
		{
			for (int f = 0; f < obj->numFaces; ++f)
			{
				if (obj->textures[f] == nullptr)
					break;

				if (!obj->faceVisible[f])
					continue;
				
				const auto [el, er] = obj->faceEdges[f];
				const auto [il, dl] = el;
				const auto [ir, dr] = er;


				std::cout << "Indx: " << il << ", " << ir << "\t\t";
				std::cout << "Dist: " << dl << ", " << dr << std::endl;
				const float xl = il * ScreenWidth() / float(frog->numRays);
				const float xr = ir * ScreenWidth() / float(frog->numRays);
				const float hl = ScreenHeight() / dl;
				const float hr = ScreenHeight() / dr;

				DrawWarpedDecal(carDecal, {
					{ xl, mY - hl / 2 },
					{ xl, mY + hl / 2 },
					{ xr, mY + hr / 2 },
					{ xr, mY - hr / 2 },
				}, obj->color);
			}
		}
		*/

		
		/*2D minimap*/
		// bg
		FillRect(0, 0,
			int(worldWidth * minimapScale),
			int(worldHeight * minimapScale), olc::BLACK);

		// rays
		for (const Player::Ray& ray : frog->rays)
			if (ray.hit != NOVD2D)
				DrawLine(
					frog->position * minimapScale,
					ray.hit * minimapScale, olc::VERY_DARK_BLUE);

		// player
		Draw(frog->position * minimapScale, olc::DARK_BLUE);

		// world objects
		for (WorldObject* obj : objs)
			for (const auto& [p1, p2] : obj->polygon)
				DrawLine(
					(obj->position + p1) * minimapScale,
					(obj->position + p2) * minimapScale, obj->color);


		// frametime
		//DrawString({ 0, ScreenHeight() - 10 }, std::to_string(fElapsedTime));

		return true;
	}
};

int main()
{
	olcWolfenFrogger game;
	if (game.Construct(200, 150, 4, 4))
		game.Start();
	return 0;
}
