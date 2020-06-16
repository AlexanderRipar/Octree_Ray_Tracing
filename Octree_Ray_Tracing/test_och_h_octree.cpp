#include "test_och_h_octree.h"

#include "och_h_octree.h"

#include <cstdint>
#include <random>
#include <chrono>
#include <cstdio>

#include "och_vec.h"
#include "och_tree_helper.h"
#include "och_string_util.h"

#include "opensimplex.h"


#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

typedef och::h_octree<19, 12> tree_t;

constexpr int screen_size_x = 640 * 2, screen_size_y = 360 * 2;
constexpr int pixel_size_x = 1, pixel_size_y = 1;

const olc::Pixel colours[8] {

	{ 0xE5, 0x8F, 0x27 },		//x_pos
	{ 0x4E, 0x61, 0x11 },		//y_pos
	{ 0x63, 0x75, 0x13 },		//z_pos
	{ 0x3D, 0x56, 0x20 },		//x_neg
	{ 0x4E, 0x61, 0x11 },		//y_neg
	{ 0x63, 0x75, 0x13 },		//z_neg
	{ 0x00, 0xBF, 0xFE },		//exit
	{ 0x3F, 0x19, 0x07 },		//inside

	//Contrast colour-scheme
	//{ 0xFF, 0x00, 0x00 },		//x_pos
	//{ 0x7F, 0x00, 0x00 },		//x_neg
	//{ 0x00, 0xFF, 0x00 },		//y_pos
	//{ 0x00, 0x7F, 0x00 },		//y_neg
	//{ 0x00, 0x00, 0xFF },		//z_pos
	//{ 0x00, 0x00, 0x7F },		//z_neg
	//{ 0x3F, 0x19, 0x07 },		//inside
	//{ 0x77, 0xFF, 0xB9 },		//exit
	//{ 0xE1, 0x07, 0xBC }		//error
};

OpenSimplexNoise terrain_noise(8789);


struct tree_camera
{
	static constexpr float delta_dir = 1.0F;

	static constexpr float delta_speed = 0.00390625F;						//1/256

	float speed = 0.0078125F;												//1/128

	const tree_t& tree;

	och::float3* rays = new och::float3[screen_size_x * screen_size_y];
	olc::Pixel* img = new olc::Pixel[screen_size_x * screen_size_y];

	bool is_flying_camera = true;

	och::float2 dir{ 0.0F, 0.0F };

	och::float3 pos{ 1.5F ,1.5F, 1.5F };

	tree_camera(const tree_t& tree) : tree{ tree } {}

	~tree_camera()
	{
		delete[] rays;
		delete[] img;
	}

	void trace_pixel(size_t x, size_t y)
	{
		size_t idx = (x + y * screen_size_x);

		och::direction hit_direction;

		uint32_t hit_voxel;

		float hit_time;


		tree.sse_trace(pos, rays[idx], hit_direction, hit_voxel, hit_time);

		//img[idx] = hit_time == 0.0F ? olc::BLACK : olc::Pixel( (uint8_t)(hit_time * 255), (uint8_t)(hit_time * 1023), (uint8_t)(hit_time * 2047));

		img[idx] = colours[static_cast<int>(hit_direction)];
	}

	void update_position()
	{
		constexpr float aspect_ratio_factor = static_cast<float>(screen_size_x) / static_cast<float>(screen_size_y);

		constexpr float view_factor_x = 2.0F / static_cast<float>(screen_size_x);

		constexpr float view_factor_y = 2.0F / static_cast<float>(screen_size_y);

		constexpr float fov = 1.25F;

		const float fov_factor = 1 / tanf(fov / 2);

		//precalculate rotation matrix-stuff
		const float sin_a = 0;
		const float cos_a = 1;
		const float sin_b = sinf(dir.x);
		const float cos_b = cosf(dir.x);
		const float sin_c = sinf(dir.y);
		const float cos_c = cosf(dir.y);

		const float t_x_fx = cos_a * cos_b;
		const float t_x_fy = cos_a * sin_b * sin_c - sin_a * cos_c;
		const float t_x_fz = cos_a * sin_b * cos_c + sin_a * sin_c;
		const float t_y_fx = sin_a * cos_b;
		const float t_y_fy = sin_a * sin_b * sin_c + cos_a * cos_c;
		const float t_y_fz = sin_a * sin_b * cos_c - cos_a * sin_c;
		const float t_z_fx = -sin_b;
		const float t_z_fy = cos_b * sin_c;
		const float t_z_fz = cos_b * cos_c;

		int idx = 0;

		for (int vert = 0; vert < screen_size_y; ++vert)
		{
			for (int horiz = 0; horiz < screen_size_x; ++horiz)
			{
				const float u = aspect_ratio_factor * (view_factor_x * static_cast<float>(horiz) - 1.0F);

				const float v = view_factor_y * static_cast<float>(vert) - 1.0F;

				const float w = fov_factor;

				const float ru = u * t_x_fx + v * t_x_fy + fov_factor * t_x_fz;
				const float rv = u * t_y_fx + v * t_y_fy + fov_factor * t_y_fz;
				const float rw = u * t_z_fx + v * t_z_fy + fov_factor * t_z_fz;

				const float reciprocal_mag = 1 / sqrtf(ru * ru + rv * rv + rw * rw);

				rays[idx++] = { rw * reciprocal_mag, ru * reciprocal_mag, -rv * reciprocal_mag };
			}
		}
	}

	olc::Pixel& at(int x, int y)
	{
		return img[x + y * screen_size_x];
	}

	void draw_crosshair()
	{
		constexpr int cx = screen_size_x / 2 - 1;
		constexpr int cy = screen_size_y / 2 - 1;

		const olc::Pixel crosshair_col(70, 30, 0);

		at(cx - 1, cy - 2) = crosshair_col;
		at(cx, cy - 2) = crosshair_col;
		at(cx + 1, cy - 2) = crosshair_col;
		at(cx + 2, cy - 2) = crosshair_col;

		at(cx - 2, cy - 1) = crosshair_col;
		at(cx - 2, cy) = crosshair_col;
		at(cx - 2, cy + 1) = crosshair_col;
		at(cx - 2, cy + 2) = crosshair_col;

		at(cx - 1, cy + 3) = crosshair_col;
		at(cx, cy + 3) = crosshair_col;
		at(cx + 1, cy + 3) = crosshair_col;
		at(cx + 2, cy + 3) = crosshair_col;

		at(cx + 3, cy - 1) = crosshair_col;
		at(cx + 3, cy) = crosshair_col;
		at(cx + 3, cy + 1) = crosshair_col;
		at(cx + 3, cy + 2) = crosshair_col;
	}

	long long update_image()
	{
		std::chrono::steady_clock::time_point beg = std::chrono::steady_clock::now();

		for (size_t y = 0; y < screen_size_y; ++y)
			for (size_t x = 0; x < screen_size_x; ++x)
			{
				trace_pixel(x, y);
			}

		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

		draw_crosshair();

		return std::chrono::duration_cast<std::chrono::milliseconds> (end - beg).count();
	}
};

class tree_window : public olc::PixelGameEngine
{
public:

	tree_camera camera;

	tree_t& tree;

	tree_window(tree_t& tree) : tree{ tree }, camera{ tree } { sAppName = "VVV - (Vx Volume Visualisation)"; }

	bool OnUserCreate() override
	{
		return true;
	}

	void update_camera_setup()
	{
		if (GetKey(olc::Key::C).bPressed)
			camera.is_flying_camera = !camera.is_flying_camera;

		if (GetKey(olc::Key::NP_ADD).bPressed)
		{
			camera.speed += camera.delta_speed;
			std::cout << "camera-movement-speed incremented to: " << camera.speed << '\n';
		}
		else if (GetKey(olc::Key::NP_SUB).bPressed)
		{
			camera.speed -= camera.delta_speed;
			std::cout << "camera-movement-speed decremented to: " << camera.speed << '\n';
		}
	}

	void update_camera_pos(float fElapsedTime, const och::float3& dir3)
	{
		float horizontal_view_dir_fct = 1 / sqrtf(dir3.x * dir3.x + dir3.y * dir3.y);

		if (camera.is_flying_camera)
		{
			if (GetKey(olc::Key::W).bHeld)
			{
				camera.pos += dir3 * camera.speed * fElapsedTime;
			}
			if (GetKey(olc::Key::S).bHeld)
			{
				camera.pos -= dir3 * camera.speed * fElapsedTime;
			}
			if (GetKey(olc::Key::A).bHeld)
			{
				camera.pos.x += dir3.y * camera.speed * fElapsedTime * horizontal_view_dir_fct;
				camera.pos.y -= dir3.x * camera.speed * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::D).bHeld)
			{
				camera.pos.x -= dir3.y * camera.speed * fElapsedTime * horizontal_view_dir_fct;
				camera.pos.y += dir3.x * camera.speed * fElapsedTime * horizontal_view_dir_fct;
			}
		}
		else
		{
			if (GetKey(olc::Key::W).bHeld)
			{
				camera.pos.x += dir3.x * camera.speed * fElapsedTime * horizontal_view_dir_fct;
				camera.pos.y += dir3.y * camera.speed * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::S).bHeld)
			{
				camera.pos.x -= dir3.x * camera.speed * fElapsedTime * horizontal_view_dir_fct;
				camera.pos.y -= dir3.y * camera.speed * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::A).bHeld)
			{
				camera.pos.x += dir3.y * camera.speed * fElapsedTime * horizontal_view_dir_fct;
				camera.pos.y -= dir3.x * camera.speed * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::D).bHeld)
			{
				camera.pos.x -= dir3.y * camera.speed * fElapsedTime * horizontal_view_dir_fct;
				camera.pos.y += dir3.x * camera.speed * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::SPACE).bHeld) camera.pos.z += camera.speed * fElapsedTime;
			if (GetKey(olc::Key::SHIFT).bHeld) camera.pos.z -= camera.speed * fElapsedTime;
		}
	}

	void update_camera_dir(float fElapsedTime)
	{
		if (GetKey(olc::Key::LEFT).bHeld)  camera.dir.x -= camera.delta_dir * fElapsedTime;
		if (GetKey(olc::Key::RIGHT).bHeld) camera.dir.x += camera.delta_dir * fElapsedTime;
		if (GetKey(olc::Key::UP).bHeld)    camera.dir.y += camera.delta_dir * fElapsedTime;
		if (GetKey(olc::Key::DOWN).bHeld)  camera.dir.y -= camera.delta_dir * fElapsedTime;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		och::float3 dir3 = { cosf(camera.dir.x) * cosf(camera.dir.y), sinf(camera.dir.x) * cosf(camera.dir.y), sinf(camera.dir.y) };

		update_camera_setup();

		update_camera_pos(fElapsedTime, dir3);

		update_camera_dir(fElapsedTime);

		float collision_dist;

		och::direction direction_dump;

		uint32_t voxel_dump;

		tree.sse_trace(camera.pos, dir3, direction_dump, voxel_dump, collision_dist);

		if (GetMouse(2).bPressed)//If the middle mouse button is pressed...
		{
			std::cout << "position (" << camera.pos.x << ", " << camera.pos.y << ", " << camera.pos.z << ")\n";
			std::cout << "facing (" << dir3.x << ", " << dir3.y << ", " << dir3.z << ")\n";

			constexpr float min_jump_dist = 0.0625F;

			collision_dist = collision_dist > min_jump_dist ? collision_dist - min_jump_dist : 0;

			camera.pos += dir3 * collision_dist;
		}
		else if (GetMouse(0).bPressed)//If left mouse button is pressed...
		{
			//constexpr float max_interact_dist = 0.0078125F;
			//
			//if (collision_dist <= max_interact_dist)
			//{
			//	och::float3 collision_pos = camera.pos + (dir3 * (collision_dist + 0.0001F));
			//	collision_pos.x -= 1.0F;
			//	collision_pos.y -= 1.0F;
			//	collision_pos.z -= 1.0F;
			//
			//	uint16_t collision_x = static_cast<uint16_t>(collision_pos.x * tree.dim);
			//	uint16_t collision_y = static_cast<uint16_t>(collision_pos.y * tree.dim);
			//	uint16_t collision_z = static_cast<uint16_t>(collision_pos.z * tree.dim);
			//
			//	tree.unset(collision_x, collision_y, collision_z);
			//}
		}
		else if (GetMouse(1).bPressed)//If right mouse button is pressed
			printf("\nPressed Right\n");

		camera.update_position();

		long long trace_time = camera.update_image();

		for (int y = 0; y < ScreenHeight(); ++y)
			for (int x = 0; x < ScreenWidth(); ++x)
				Draw(x, y, olc::Pixel(camera.img[x + y * ScreenWidth()]));

		DrawString(8, 8, std::to_string(trace_time) + " ms");

		return true;
	}
};



int get_terrain_heigth(int x, int y, size_t dim)
{
	return static_cast<int>(terrain_noise.Evaluate((static_cast<float>(x * 4) / dim), (static_cast<float>(y * 4) / dim)) * dim / 16 + dim / 2);
}

struct heightmap
{
	size_t dim;

	uint16_t* data = new uint16_t[dim * dim];

	uint16_t at(int x, int y) const
	{
		return data[y * dim + x];
	}

	uint16_t& get(int x, int y)
	{
		return data[y * dim + x];
	}

	heightmap(int dim) : dim(dim)
	{
		for(int y = 0; y < dim; ++y)
			for (int x = 0; x < dim; ++x)
				get(x, y) = get_terrain_heigth(x, y, dim);
	}

	~heightmap()
	{
		delete[] data;
	}
};

uint32_t get_leaf_val()
{
	return 1; //(std::rand() & 0b1111111111) + 1;
}

uint32_t create_node(const int x, const int y, const int z, const int depth, const heightmap& heights, tree_t& tree)
{
	const int dim = 1 << depth;

	tree_t::node node;

	for(int _y = 0; _y < dim; ++_y)
		for (int _x = 0; _x < dim; ++_x)
		{
			uint16_t height = heights.at(x + _x, y + _y);

			if (height >= z && height < z + dim)
				goto ACTIVE;
		}

	return 0;

ACTIVE:

	const int half_dim = 1 << (depth - 1);

	if (depth != 1)//No leaf yet...
	{
		node.children[0] = create_node(x           , y           , z           , depth - 1, heights, tree);
		node.children[1] = create_node(x + half_dim, y           , z           , depth - 1, heights, tree);
		node.children[2] = create_node(x           , y + half_dim, z           , depth - 1, heights, tree);
		node.children[3] = create_node(x + half_dim, y + half_dim, z           , depth - 1, heights, tree);
		node.children[4] = create_node(x           , y           , z + half_dim, depth - 1, heights, tree);
		node.children[5] = create_node(x + half_dim, y           , z + half_dim, depth - 1, heights, tree);
		node.children[6] = create_node(x           , y + half_dim, z + half_dim, depth - 1, heights, tree);
		node.children[7] = create_node(x + half_dim, y + half_dim, z + half_dim, depth - 1, heights, tree);
	}
	else//Leaf
	{
		node.children[0] = heights.at(x    , y    ) == z     ? get_leaf_val() : 0;
		node.children[1] = heights.at(x + 1, y    ) == z     ? get_leaf_val() : 0;
		node.children[2] = heights.at(x    , y + 1) == z     ? get_leaf_val() : 0;
		node.children[3] = heights.at(x + 1, y + 1) == z     ? get_leaf_val() : 0;
		node.children[4] = heights.at(x    , y    ) == z + 1 ? get_leaf_val() : 0;
		node.children[5] = heights.at(x + 1, y    ) == z + 1 ? get_leaf_val() : 0;
		node.children[6] = heights.at(x    , y + 1) == z + 1 ? get_leaf_val() : 0;
		node.children[7] = heights.at(x + 1, y + 1) == z + 1 ? get_leaf_val() : 0;
	}

	if (node.is_zero())
		printf("ERR");

	return tree.register_node(node);
}



void test_och_h_octree()
{
	printf("\nEntered test_och_h_octree...\n");

	tree_t tree;

	tree_window window(tree);

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	{
		heightmap heightfield(tree.dim);
	
		tree.set_root(create_node(0, 0, 0, tree.depth, heightfield, tree));
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	printf("\nTree-depth: %i\n", tree.depth);

	printf("\nTree-dimension: %i\n", tree.dim);

	printf("\nTime taken to construct svodag:   %lli ms\n", std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count());

	printf("\nFilled entries in node-hashtable: %i of %i max\n", tree.get_fillcnt(), tree.table_capacity);

	printf("\nTotal node-count: %i\n", tree.get_nodecnt());

	printf("\nCompression-ratio: %f\n", (float)tree.get_nodecnt() / tree.get_fillcnt());

	printf("\nRequired size in memory: %s\n", och::abbreviate_byte_size(static_cast<size_t>(tree.get_fillcnt()) * 37).c_str());

	printf("\nUncompressed size: %s\n", och::abbreviate_byte_size(static_cast<size_t>(tree.get_nodecnt()) * 32).c_str());

	printf("\nTable-size: %s\n", och::abbreviate_byte_size(tree.table_bytes).c_str());

	if (window.Construct(screen_size_x, screen_size_y, pixel_size_x, pixel_size_y))
		window.Start();

	printf("\nExited test_och_hashed_octree\n");
}
