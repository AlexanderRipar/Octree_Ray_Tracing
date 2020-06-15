#include "test_och_h_octree.h"

#include <cstdint>
#include <random>
#include <chrono>
#include <random>

#include "och_vec.h"
#include "och_tree_helper.h"
#include "och_h_octree.h"
#include "och_string_util.h"

#include "opensimplex.h"


#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

constexpr bool is_norm = true;

constexpr int screen_size_x = 640, screen_size_y = 360;
constexpr int pixel_size_x = 2, pixel_size_y = 2;

constexpr int octree_depth = 12;
constexpr int octree_dim = 1 << octree_depth;

constexpr olc::Pixel colours[8] {

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

constexpr olc::Pixel voxel_colours[2] {
	{ 0x77, 0xFF, 0xB9 },		//air
	{ 0x3F, 0x19, 0x07 }		//active
};

OpenSimplexNoise terrain_noise(8789);





struct h_octree_camera
{
	och::float3 rays[screen_size_x * screen_size_y];
	olc::Pixel img[screen_size_x * screen_size_y];

	och::float3 camera_pos{ 1,1,1 };

	void trace_pixel(size_t x, size_t y)
	{
		size_t idx = (x + y * screen_size_x);

		och::direction hit_direction;

		och::voxel hit_voxel;

		float hit_time;


		if constexpr(is_norm)
			och::sse_trace_h_octree(octree_depth, camera_pos, rays[idx], hit_direction, hit_voxel, hit_time);
		else
		{
			och::inv_trace_h_octree(octree_depth - 1, camera_pos, rays[idx], hit_direction, hit_voxel, hit_time);
		}

		img[idx] = colours[static_cast<int>(hit_direction)];
	}

	void update_position(och::float2 v_angles, och::float3 v_pos)
	{
		constexpr float aspect_ratio_factor = static_cast<float>(screen_size_x) / static_cast<float>(screen_size_y);

		constexpr float view_factor_x = 2.0F / static_cast<float>(screen_size_x);

		constexpr float view_factor_y = 2.0F / static_cast<float>(screen_size_y);

		constexpr float fov = 1.25F;

		const float fov_factor = 1 / tanf(fov / 2);

		//precalculate rotation matrix-stuff
		const float sin_a = 0;
		const float cos_a = 1;
		const float sin_b = sinf(v_angles.x);
		const float cos_b = cosf(v_angles.x);
		const float sin_c = sinf(v_angles.y);
		const float cos_c = cosf(v_angles.y);

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

		if constexpr(is_norm)
		{
			camera_pos.x = v_pos.x / octree_dim + 1.0F;
			camera_pos.y = v_pos.y / octree_dim + 1.0F;
			camera_pos.z = v_pos.z / octree_dim + 1.0F;
		}
		else
			camera_pos = v_pos;
	}

	void draw_crosshair()
	{
		int center_x = screen_size_x / 2 - 1;
		int center_y = screen_size_y / 2 - 1;

		int center = center_x + center_y * screen_size_x;

		img[center - 5] = { 0, 0, 0 };
		img[center - 4] = { 0, 0, 0 };
		img[center - 3] = { 0, 0, 0 };
		img[center - 2] = { 0, 0, 0 };
		img[center - 1] = { 0, 0, 0 };
		img[center] = { 0, 0, 0 };
		img[center + 1] = { 0, 0, 0 };
		img[center + 2] = { 0, 0, 0 };
		img[center + 3] = { 0, 0, 0 };
		img[center + 4] = { 0, 0, 0 };
		img[center + 5] = { 0, 0, 0 };
		img[center + 6] = { 0, 0, 0 };

		img[center + screen_size_x - 5] = { 0, 0, 0 };
		img[center + screen_size_x - 4] = { 0, 0, 0 };
		img[center + screen_size_x - 3] = { 0, 0, 0 };
		img[center + screen_size_x - 2] = { 0, 0, 0 };
		img[center + screen_size_x - 1] = { 0, 0, 0 };
		img[center + screen_size_x] = { 0, 0, 0 };
		img[center + screen_size_x + 1] = { 0, 0, 0 };
		img[center + screen_size_x + 2] = { 0, 0, 0 };
		img[center + screen_size_x + 3] = { 0, 0, 0 };
		img[center + screen_size_x + 4] = { 0, 0, 0 };
		img[center + screen_size_x + 5] = { 0, 0, 0 };
		img[center + screen_size_x + 6] = { 0, 0, 0 };

		img[center - screen_size_x * 5] = { 0, 0, 0 };
		img[center - screen_size_x * 4] = { 0, 0, 0 };
		img[center - screen_size_x * 3] = { 0, 0, 0 };
		img[center - screen_size_x * 2] = { 0, 0, 0 };
		img[center - screen_size_x] = { 0, 0, 0 };
		img[center + screen_size_x * 2] = { 0, 0, 0 };
		img[center + screen_size_x * 3] = { 0, 0, 0 };
		img[center + screen_size_x * 4] = { 0, 0, 0 };
		img[center + screen_size_x * 5] = { 0, 0, 0 };
		img[center + screen_size_x * 6] = { 0, 0, 0 };

		img[center - screen_size_x * 5 + 1] = { 0, 0, 0 };
		img[center - screen_size_x * 4 + 1] = { 0, 0, 0 };
		img[center - screen_size_x * 3 + 1] = { 0, 0, 0 };
		img[center - screen_size_x * 2 + 1] = { 0, 0, 0 };
		img[center - screen_size_x + 1] = { 0, 0, 0 };
		img[center + screen_size_x * 2 + 1] = { 0, 0, 0 };
		img[center + screen_size_x * 3 + 1] = { 0, 0, 0 };
		img[center + screen_size_x * 4 + 1] = { 0, 0, 0 };
		img[center + screen_size_x * 5 + 1] = { 0, 0, 0 };
		img[center + screen_size_x * 6 + 1] = { 0, 0, 0 };


		img[center + screen_size_x + 1] = { 0, 0, 0 };
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

h_octree_camera camera;

class octnode_window_test : public olc::PixelGameEngine
{
public:

	float camera_view_pos_change = 16.0F;

	och::float2 camera_view_dir{ 0, 0 };

	och::float3 camera_view_pos{ 1 << (octree_depth - 1), 1 << (octree_depth - 1) , 1 << (octree_depth - 1) };

	bool is_directional_flying_camera = true;

	octnode_window_test() { sAppName = "VVV - (Vx Volume Visualisation)"; }

	bool OnUserCreate() override
	{
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		constexpr float view_pos_incr_change = 1.0F;

		if (GetKey(olc::Key::C).bPressed) is_directional_flying_camera = !is_directional_flying_camera;

		if (GetKey(olc::Key::NP_ADD).bPressed)
		{
			camera_view_pos_change += view_pos_incr_change;
			std::cout << "camera-movement-speed incremented to: " << camera_view_pos_change << '\n';
		}
		else if (GetKey(olc::Key::NP_SUB).bPressed)
		{
			camera_view_pos_change -= view_pos_incr_change;
			std::cout << "camera-movement-speed decremented to: " << camera_view_pos_change << '\n';
		}

		camera.update_position(camera_view_dir, camera_view_pos);

		och::float3 center_ray = camera.rays[screen_size_x / 2 + (screen_size_y / 2) * screen_size_x];

		float horizontal_view_dir_fct = 1 / sqrtf(center_ray.x * center_ray.x + center_ray.y * center_ray.y);

		if (is_directional_flying_camera)
		{
			if (GetKey(olc::Key::W).bHeld)
			{
				camera_view_pos.x += center_ray.x * camera_view_pos_change * fElapsedTime;
				camera_view_pos.y += center_ray.y * camera_view_pos_change * fElapsedTime;
				camera_view_pos.z += center_ray.z * camera_view_pos_change * fElapsedTime;
			}
			if (GetKey(olc::Key::S).bHeld)
			{
				camera_view_pos.x -= center_ray.x * camera_view_pos_change * fElapsedTime;
				camera_view_pos.y -= center_ray.y * camera_view_pos_change * fElapsedTime;
				camera_view_pos.z -= center_ray.z * camera_view_pos_change * fElapsedTime;
			}
			if (GetKey(olc::Key::A).bHeld)
			{
				camera_view_pos.x += center_ray.y * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
				camera_view_pos.y -= center_ray.x * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::D).bHeld)
			{
				camera_view_pos.x -= center_ray.y * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
				camera_view_pos.y += center_ray.x * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
			}
		}
		else
		{
			if (GetKey(olc::Key::W).bHeld)
			{
				camera_view_pos.x += center_ray.x * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
				camera_view_pos.y += center_ray.y * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::S).bHeld)
			{
				camera_view_pos.x -= center_ray.x * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
				camera_view_pos.y -= center_ray.y * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::A).bHeld)
			{
				camera_view_pos.x += center_ray.y * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
				camera_view_pos.y -= center_ray.x * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::D).bHeld)
			{
				camera_view_pos.x -= center_ray.y * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
				camera_view_pos.y += center_ray.x * camera_view_pos_change * fElapsedTime * horizontal_view_dir_fct;
			}
			if (GetKey(olc::Key::SPACE).bHeld) camera_view_pos.z += camera_view_pos_change * fElapsedTime;
			if (GetKey(olc::Key::SHIFT).bHeld) camera_view_pos.z -= camera_view_pos_change * fElapsedTime;
		}

		constexpr float view_direction_speed = 1.0F;

		if (GetKey(olc::Key::LEFT).bHeld) camera_view_dir.x -= view_direction_speed * fElapsedTime;
		if (GetKey(olc::Key::RIGHT).bHeld) camera_view_dir.x += view_direction_speed * fElapsedTime;
		if (GetKey(olc::Key::UP).bHeld) camera_view_dir.y += view_direction_speed * fElapsedTime;
		if (GetKey(olc::Key::DOWN).bHeld) camera_view_dir.y -= view_direction_speed * fElapsedTime;

		if (GetMouse(0).bPressed)//If the left mouse button is pressed...
		{
			float view_dir_collision_dist;

			och::float3 view_dir_3 = { cosf(camera_view_dir.x) * cosf(camera_view_dir.y), sinf(camera_view_dir.x) * cosf(camera_view_dir.y), sinf(camera_view_dir.y) };

			std::cout << "position (" << camera_view_pos.x << ", " << camera_view_pos.y << ", " << camera_view_pos.z << ")\n";
			std::cout << "facing (" << view_dir_3.x << ", " << view_dir_3.y << ", " << view_dir_3.z << ")\n";

			och::direction dir_dump;

			och::voxel vox_dump;

			if constexpr (is_norm)
				och::sse_trace_h_octree(octree_depth, camera_view_pos, view_dir_3, dir_dump, vox_dump, view_dir_collision_dist);
			else
				och::inv_trace_h_octree(octree_depth - 1, camera_view_pos, view_dir_3, dir_dump, vox_dump, view_dir_collision_dist);

			view_dir_collision_dist = view_dir_collision_dist > 32 ? view_dir_collision_dist - 32 : 0;

			camera_view_pos += view_dir_3 * view_dir_collision_dist;
		}

		long long trace_time = camera.update_image();

		for (int y = 0; y < ScreenHeight(); ++y)
			for (int x = 0; x < ScreenWidth(); ++x)
				Draw(x, y, olc::Pixel(camera.img[x + y * ScreenWidth()]));

		DrawString(8, 8, std::to_string(trace_time) + " ms");

		return true;
	}
};

octnode_window_test window;





int get_terrain_heigth(int x, int y)
{
	return static_cast<int>(terrain_noise.Evaluate((static_cast<float>(x * 4) / octree_dim), (static_cast<float>(y * 4) / octree_dim)) * octree_dim / 16 + octree_dim / 2);
}

uint32_t create_node_try4_inline(int x, int y, int z, int depth)
{
	och::h_node_8 node;

	int nx = x << 1, ny = y << 1, nz = z << 1;

	if (depth == 1)
	{
		node.children[0] = get_terrain_heigth(nx    , ny    ) == nz    ;
		node.children[1] = get_terrain_heigth(nx + 1, ny    ) == nz    ;
		node.children[2] = get_terrain_heigth(nx    , ny + 1) == nz    ;
		node.children[3] = get_terrain_heigth(nx + 1, ny + 1) == nz    ;
		node.children[4] = get_terrain_heigth(nx    , ny    ) == nz + 1;
		node.children[5] = get_terrain_heigth(nx + 1, ny    ) == nz + 1;
		node.children[6] = get_terrain_heigth(nx    , ny + 1) == nz + 1;
		node.children[7] = get_terrain_heigth(nx + 1, ny + 1) == nz + 1;
	}
	else
	{
		node.children[0] = create_node_try4_inline(nx    , ny    , nz    , depth - 1);
		node.children[1] = create_node_try4_inline(nx + 1, ny    , nz    , depth - 1);
		node.children[2] = create_node_try4_inline(nx    , ny + 1, nz    , depth - 1);
		node.children[3] = create_node_try4_inline(nx + 1, ny + 1, nz    , depth - 1);
		node.children[4] = create_node_try4_inline(nx    , ny    , nz + 1, depth - 1);
		node.children[5] = create_node_try4_inline(nx + 1, ny    , nz + 1, depth - 1);
		node.children[6] = create_node_try4_inline(nx    , ny + 1, nz + 1, depth - 1);
		node.children[7] = create_node_try4_inline(nx + 1, ny + 1, nz + 1, depth - 1);
	}

	if (node == och::h_node_8{ 0, 0, 0, 0, 0, 0, 0, 0 })
		return 0;

	return register_h_node_8(&node);
}

struct tree_initializer
{
	uint16_t* data = new uint16_t[octree_dim * octree_dim];

	uint16_t at(int x, int y) const
	{
		return data[y * octree_dim + x];
	}

	uint16_t& get(int x, int y)
	{
		return data[y * octree_dim + x];
	}

	tree_initializer()
	{
		for(int y = 0; y < octree_dim; ++y)
			for (int x = 0; x < octree_dim; ++x)
				get(x, y) = get_terrain_heigth(x, y);
	}

	~tree_initializer()
	{
		delete[] data;
	}
};

uint32_t get_leaf_val()
{
	return 1; //(std::rand() & 0b1111111111) + 1;
}

uint32_t create_node_try7(const int x, const int y, const int z, const int depth, const tree_initializer& heights)
{
	const int dim = 1 << depth;

	och::h_node_8 node;

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
		node.children[0] = create_node_try7(x           , y           , z           , depth - 1, heights);
		node.children[1] = create_node_try7(x + half_dim, y           , z           , depth - 1, heights);
		node.children[2] = create_node_try7(x           , y + half_dim, z           , depth - 1, heights);
		node.children[3] = create_node_try7(x + half_dim, y + half_dim, z           , depth - 1, heights);
		node.children[4] = create_node_try7(x           , y           , z + half_dim, depth - 1, heights);
		node.children[5] = create_node_try7(x + half_dim, y           , z + half_dim, depth - 1, heights);
		node.children[6] = create_node_try7(x           , y + half_dim, z + half_dim, depth - 1, heights);
		node.children[7] = create_node_try7(x + half_dim, y + half_dim, z + half_dim, depth - 1, heights);
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

	och::register_h_node_8(&node);
}

void test_och_hashed_octree()
{
	std::cout << "test_och_hashed_octree...\n";

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	{
		tree_initializer heightfield;
	
		och::set_root(create_node_try7(0, 0, 0, octree_depth, heightfield));
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	printf("\nTree-depth: %i\n", octree_depth);

	printf("\nTree-dimension: %i\n", octree_dim);

	printf("\nTime taken to construct svodag:   %lli ms\n", std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count());

	printf("\nFilled entries in node-hashtable: %i of %i max\n", och::get_fill_count(), och::get_table_size());

	printf("\nTotal node-count: %i\n", och::get_node_count());

	printf("\nCompression-ratio: %f (%i%%)\n", (float)och::get_node_count() / och::get_fill_count(), (int)(100.0F / ((float)och::get_fill_count() / och::get_node_count())));

	printf("\nRequired size in memory: %s\n", och::abbreviate_byte_size(static_cast<size_t>(och::get_fill_count()) * 35).c_str());

	printf("\nUncompressed size: %s\n", och::abbreviate_byte_size(static_cast<size_t>(och::get_node_count()) * 32).c_str());

	if (window.Construct(screen_size_x, screen_size_y, pixel_size_x, pixel_size_y))
		window.Start();

	std::cout << "\nexited test_och_hashed_octree\n";
}
