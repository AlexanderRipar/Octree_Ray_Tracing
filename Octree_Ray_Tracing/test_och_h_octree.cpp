#include "test_och_h_octree.h"

#include "och_h_octree.h"

#include <cstdint>
#include <random>
#include <chrono>
#include <cstdio>

#include <Windows.h>

#include "och_vec.h"
#include "och_tree_helper.h"
#include "och_string_util.h"
#include "och_voxel.h"

#include "och_noise.h"
#include "opensimplex.h"


#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

typedef och::h_octree<19, 8> tree_t;

constexpr bool is_half_res = true;

constexpr int screen_size_x = 640 * (2 - is_half_res);
constexpr int screen_size_y = 360 * (2 - is_half_res);
constexpr int pixel_size  = 1 + is_half_res;

OpenSimplexNoise terrain_noise(8789);

och::simplex_n noise(0.5F);

struct tree_camera
{
	static constexpr float delta_dir_mouse = 1.0F / 48;

	static constexpr float delta_speed = 1.0F / tree_t::dim;								//corresponds to one voxel

	float speed = delta_speed * (tree_t::dim / 32);

	const tree_t& tree;

	const och::voxel_data& voxels;

	och::float3* rays = new och::float3[screen_size_x * screen_size_y];

	bool is_flying_camera = true;

	och::float2 dir{ 0.0F, 0.0F };

	och::float3 pos{ 1.5F ,1.5F, 1.5F };

	tree_camera(const tree_t& tree, const och::voxel_data& voxels) : tree{ tree }, voxels{ voxels } {}

	~tree_camera()
	{
		delete[] rays;
	}

	olc::Pixel trace_pixel(size_t x, size_t y)
	{
		size_t idx = (x + y * screen_size_x);

		och::direction hit_direction;

		uint32_t hit_voxel;

		float hit_time;

		tree.sse_trace(pos, rays[idx], hit_direction, hit_voxel, hit_time);

		const olc::Pixel exit_colour{ 0x00, 0xBF, 0xFE };
		const olc::Pixel inside_colour{ 0x3F, 0x19, 0x07 };

		if (hit_direction == och::direction::exit)
			return exit_colour;
		else if (hit_direction == och::direction::inside)
			return inside_colour;

		return reinterpret_cast<const olc::Pixel*>(voxels.colours)[6 * (hit_voxel - 1) + static_cast<uint32_t>(hit_direction)];
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
};

class tree_window : public olc::PixelGameEngine
{
public:

	static constexpr int txt_beg_x = 8, txt_beg_y = 8, txt_incr = 8;

	static constexpr float max_interact_dist = 0.5F;

	tree_t& tree;

	const och::voxel_data& voxels;

	tree_camera camera;

	bool is_debug = true;

	och::vec3f measure_pt{ -1.0F, -1.0F, -1.0F };
	std::string measure_output;

	tree_window(tree_t& tree, const och::voxel_data& voxels) : tree{ tree }, voxels{ voxels }, camera{ tree, voxels } { sAppName = "VVV - (Vx Volume Visualisation)"; }

	bool OnUserCreate() override
	{
		return true;
	}

	void update_camera_setup()
	{
		if (GetKey(olc::Key::C).bPressed)
			camera.is_flying_camera = !camera.is_flying_camera;

		if(GetKey(olc::Key::SHIFT).bHeld)
			camera.speed += camera.delta_speed * (GetMouseWheel() / 120);
		else
			camera.speed += camera.delta_speed * (GetMouseWheel() / 15);

		if (camera.speed < 0)
			camera.speed = 0;
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
		RECT window_rect;

		GetWindowRect(olc_hWnd, &window_rect);

		POINT mouse_pos;

		GetCursorPos(&mouse_pos);

		SetCursorPos(window_rect.right / 2, window_rect.bottom / 2);

		int mouse_delta_x = mouse_pos.x - window_rect.right / 2;
		int mouse_delta_y = mouse_pos.y - window_rect.bottom / 2;

		camera.dir.x += mouse_delta_x * camera.delta_dir_mouse * fElapsedTime;
		camera.dir.y -= mouse_delta_y * camera.delta_dir_mouse * fElapsedTime;
	}

	void update_text(long long trace_time, const och::float3& dir3, const och::float3& offset, float hit_dst, uint32_t hit_vox, och::direction hit_dir)
	{
		if (GetKey(olc::Key::O).bPressed)
			is_debug = !is_debug;

		if (is_debug)
		{
			//Handle "Looking at"
			och::float3 collision_pos = camera.pos + (dir3 * hit_dst) + offset;

			collision_pos.x -= 1.0F;
			collision_pos.y -= 1.0F;
			collision_pos.z -= 1.0F;

			uint16_t collision_x = static_cast<uint16_t>(collision_pos.x * tree.dim);
			uint16_t collision_y = static_cast<uint16_t>(collision_pos.y * tree.dim);
			uint16_t collision_z = static_cast<uint16_t>(collision_pos.z * tree.dim);

			std::string looking_at_str;

			looking_at_str = hit_vox ? "[" + std::to_string(collision_x) + ", " + std::to_string(collision_y) + ", " + std::to_string(collision_z) + "]: " + voxels.names[hit_vox - 1].str + " (" + std::to_string(hit_vox) + ")" : "Air (0)";

			//Handle "facing"
			float dir_x = abs(dir3.x);
			float dir_y = abs(dir3.y);
			float dir_z = abs(dir3.z);

			std::string facing_str = std::to_string(dir3.x) + ", " + std::to_string(dir3.y) + ", " + std::to_string(dir3.z);

			if (dir_x >= dir_y && dir_x >= dir_z)
				facing_str += dir3.x < 0 ? " (-x)" : " (+x)";
			else if (dir_y > dir_x && dir_y >= dir_z)
				facing_str += dir3.y < 0 ? " (-y)" : " (+y)";
			else
				facing_str += dir3.z < 0 ? " (-z)" : " (+z)";

			//Output
			DrawString(txt_beg_x, txt_beg_y + txt_incr * 0, std::to_string(trace_time) + " ms");
			DrawString(txt_beg_x, txt_beg_y + txt_incr * 1, "tabled nodes: " + std::to_string(tree.get_fillcnt()));
			DrawString(txt_beg_x, txt_beg_y + txt_incr * 2, "active nodes: " + std::to_string(tree.get_nodecnt()));
			DrawString(txt_beg_x, txt_beg_y + txt_incr * 3, "memory usage: " + och::abbreviate_byte_size(static_cast<size_t>(tree.get_fillcnt()) * 37));
			DrawString(txt_beg_x, txt_beg_y + txt_incr * 4, "speed:        " + std::to_string(static_cast<int>(camera.speed * tree.dim)));
			DrawString(txt_beg_x, txt_beg_y + txt_incr * 5, "facing:       " + facing_str);
			DrawString(txt_beg_x, txt_beg_y + txt_incr * 6, "Looking at:   " + looking_at_str);
		}

		DrawString(txt_beg_x + txt_incr / 2, txt_beg_y + txt_incr * 7 + txt_incr / 2, measure_output);
	}

	void user_interaction(const och::float3& dir3, const och::float3& offset, float hit_dst)
	{
		if (GetMouse(2).bPressed)//If the middle mouse button is pressed...
		{
			std::cout << "position (" << camera.pos.x << ", " << camera.pos.y << ", " << camera.pos.z << ")\n";
			std::cout << "facing (" << dir3.x << ", " << dir3.y << ", " << dir3.z << ")\n";

			constexpr float min_jump_dist = 0.0625F;

			hit_dst = hit_dst > min_jump_dist ? hit_dst - min_jump_dist : 0;

			camera.pos += dir3 * hit_dst;
		}

		if (GetKey(olc::Key::I).bPressed)
		{
			int x = static_cast<int>((camera.pos.x - 1.0F) * tree.dim);
			int y = static_cast<int>((camera.pos.y - 1.0F) * tree.dim);
			int z = static_cast<int>((camera.pos.z - 1.0F) * tree.dim);

			int init_z = z;

			while (z != tree.dim && tree.at(x, y, z))
				z++;

			if (z != init_z && z != tree.dim)
				camera.pos.z = static_cast<float>(z + 1) / tree.dim + 1.0F;
		}

		if (GetKey(olc::Key::M).bPressed)
		{
			if (measure_output != "" || measure_pt.x == -1.0F)				//measure_0
			{
				measure_pt.x = ((camera.pos.x - 1.0F) + dir3.x * hit_dst) * tree.dim;
				measure_pt.y = ((camera.pos.y - 1.0F) + dir3.y * hit_dst) * tree.dim;
				measure_pt.z = ((camera.pos.z - 1.0F) + dir3.z * hit_dst) * tree.dim;

				measure_output = "";
			}
			else															//measure_1
			{
				och::vec3f measure(((camera.pos.x - 1.0F) + dir3.x * hit_dst) * tree.dim, ((camera.pos.y - 1.0F) + dir3.y * hit_dst) * tree.dim, ((camera.pos.z - 1.0F) + dir3.z * hit_dst) * tree.dim);

				float dist = sqrt((measure_pt.x - measure.x) * (measure_pt.x - measure.x) + (measure_pt.y - measure.y) * (measure_pt.y - measure.y) + (measure_pt.z - measure.z) * (measure_pt.z - measure.z));

				measure_output = "Dist = " + std::to_string(dist);
			}
		}

		if (hit_dst >= max_interact_dist)
			return;

		if (GetMouse(0).bPressed)//If left mouse button is pressed...
		{
			och::float3 collision_pos = camera.pos + (dir3 * hit_dst) + offset;

			collision_pos.x -= 1.0F;
			collision_pos.y -= 1.0F;
			collision_pos.z -= 1.0F;

			uint16_t collision_x = static_cast<uint16_t>(collision_pos.x * tree.dim);
			uint16_t collision_y = static_cast<uint16_t>(collision_pos.y * tree.dim);
			uint16_t collision_z = static_cast<uint16_t>(collision_pos.z * tree.dim);

			tree.set(collision_x, collision_y, collision_z, 0);
		}
		else if (GetMouse(1).bPressed)//If right mouse button is pressed
		{
			och::float3 collision_pos = camera.pos + (dir3 * hit_dst) - offset;
			collision_pos.x -= 1.0F;
			collision_pos.y -= 1.0F;
			collision_pos.z -= 1.0F;

			uint16_t collision_x = static_cast<uint16_t>(collision_pos.x * tree.dim);
			uint16_t collision_y = static_cast<uint16_t>(collision_pos.y * tree.dim);
			uint16_t collision_z = static_cast<uint16_t>(collision_pos.z * tree.dim);

			tree.set(collision_x, collision_y, collision_z, 1);
		}

		if (GetKey(olc::Key::T).bHeld)
		{
			och::float3 collision_pos = camera.pos + (dir3 * hit_dst) + offset;
			collision_pos.x -= 1.0F;
			collision_pos.y -= 1.0F;
			collision_pos.z -= 1.0F;

			uint16_t coll_x = static_cast<uint16_t>(collision_pos.x * tree.dim);
			uint16_t coll_y = static_cast<uint16_t>(collision_pos.y * tree.dim);
			uint16_t coll_z = static_cast<uint16_t>(collision_pos.z * tree.dim);

			constexpr int dim = 40;

			for (int z = -dim / 2; z < (dim + 1) / 2; ++z)
				for (int y = -dim / 2; y < (dim + 1) / 2; ++y)
					for (int x = -dim / 2; x < (dim + 1) / 2; ++x)
						tree.set(coll_x + x, coll_y + y, coll_z + z, 1);
		}

		if (GetKey(olc::Key::Z).bHeld)
		{
				och::float3 collision_pos = camera.pos + (dir3 * hit_dst) - offset;
				collision_pos.x -= 1.0F;
				collision_pos.y -= 1.0F;
				collision_pos.z -= 1.0F;

				uint16_t coll_x = static_cast<uint16_t>(collision_pos.x * tree.dim);
				uint16_t coll_y = static_cast<uint16_t>(collision_pos.y * tree.dim);
				uint16_t coll_z = static_cast<uint16_t>(collision_pos.z * tree.dim);

				constexpr int dim = 40;

				for (int z = -dim / 2; z < (dim + 1) / 2; ++z)
					for (int y = -dim / 2; y < (dim + 1) / 2; ++y)
						for (int x = -dim / 2; x < (dim + 1) / 2; ++x)
							tree.set(coll_x + x, coll_y + y, coll_z + z, 0);
		}

	}

	long long update_image()
	{
		const olc::Pixel exit_colour{ 0x00, 0xBF, 0xFE };

		std::chrono::steady_clock::time_point beg = std::chrono::steady_clock::now();

		if (!tree.get_root())
			for (int y = 0; y < ScreenHeight(); ++y)
				for (int x = 0; x < ScreenWidth(); ++x)
					Draw(x, y, exit_colour);
		else
			for (int y = 0; y < ScreenHeight(); ++y)
				for (int x = 0; x < ScreenWidth(); ++x)
					Draw(x, y, camera.trace_pixel(x, y));

		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

		draw_crosshair();

		return std::chrono::duration_cast<std::chrono::milliseconds> (end - beg).count();
	}

	void draw_crosshair()
	{
		constexpr int cx = screen_size_x / 2 - 1;
		constexpr int cy = screen_size_y / 2 - 1;

		const olc::Pixel crosshair_col(70, 30, 0);

		Draw(cx - 1, cy - 2, crosshair_col);
		Draw(cx    , cy - 2, crosshair_col);
		Draw(cx + 1, cy - 2, crosshair_col);
		Draw(cx + 2, cy - 2, crosshair_col);

		Draw(cx - 2, cy - 1, crosshair_col);
		Draw(cx - 2, cy    , crosshair_col);
		Draw(cx - 2, cy + 1, crosshair_col);
		Draw(cx - 2, cy + 2, crosshair_col);

		Draw(cx - 1, cy + 3, crosshair_col);
		Draw(cx    , cy + 3, crosshair_col);
		Draw(cx + 1, cy + 3, crosshair_col);
		Draw(cx + 2, cy + 3, crosshair_col);

		Draw(cx + 3, cy - 1, crosshair_col);
		Draw(cx + 3, cy    , crosshair_col);
		Draw(cx + 3, cy + 1, crosshair_col);
		Draw(cx + 3, cy + 2, crosshair_col);
	}

	och::float3 get_directional_hit_offset(och::direction hit_dir)
	{
		och::float3 offset(0, 0, 0);

		switch (hit_dir)
		{
		case och::direction::x_pos: offset.x = tree.voxel_dim / 2; break;
		case och::direction::y_pos: offset.y = tree.voxel_dim / 2; break;
		case och::direction::z_pos: offset.z = tree.voxel_dim / 2; break;
		case och::direction::x_neg: offset.x = -tree.voxel_dim / 2; break;
		case och::direction::y_neg: offset.y = -tree.voxel_dim / 2; break;
		case och::direction::z_neg: offset.z = -tree.voxel_dim / 2; break;
		}

		return offset;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		//for (int y = 0; y < ScreenHeight(); ++y)
		//	for (int x = 0; x < ScreenWidth(); ++x)
		//	{
		//		float n = noise(x, y, 4) + 1.0F;
		//
		//		Draw(x, y, olc::Pixel(0, 0, (uint8_t)(128 * n)));
		//	}

		//for (int y = 0; y < ScreenHeight(); ++y)
		//	for (int x = 0; x < ScreenWidth(); ++x)
		//	{
		//		float n = terrain_noise.Evaluate(x / 32.0F, y / 32.0F, 4) + 1.0F;
		//
		//		Draw(x, y, olc::Pixel(0, 0, (uint8_t)(128 * n)));
		//	}

		if (IsFocused())
		{
			if (GetKey(olc::Key::ESCAPE).bHeld || GetKey(olc::Key::ENTER).bHeld)
				return false;
		
			och::float3 dir3 = { cosf(camera.dir.x) * cosf(camera.dir.y), sinf(camera.dir.x) * cosf(camera.dir.y), sinf(camera.dir.y) };
		
			float hit_dst = INFINITY;
		
			och::direction hit_dir = och::direction::exit;
		
			uint32_t hit_vox = 0;
			
			if(tree.get_root())
				tree.sse_trace(camera.pos, dir3, hit_dir, hit_vox, hit_dst);
		
			och::float3 offset = get_directional_hit_offset(hit_dir);
			
			update_camera_setup();
			
			update_camera_pos(fElapsedTime, dir3);
		
			update_camera_dir(fElapsedTime);
			
			user_interaction(dir3, offset, hit_dst);
			
			camera.update_position();
			
			long long trace_time = update_image();
			
			update_text(trace_time, dir3, offset, hit_dst, hit_vox, hit_dir);
		}
		
		return true;
	}
};



int get_terrain_heigth(int x, int y)
{
	float px = (static_cast<float>(x * 4) / tree_t::dim);
	float py = (static_cast<float>(y * 4) / tree_t::dim);

	return static_cast<int>(noise(px, py) * tree_t::dim / 16 + tree_t::dim / 4);

	//return static_cast<int>(terrain_noise.Evaluate(px, py) * tree_t::dim / 16 + tree_t::dim / 4);
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
				get(x, y) = get_terrain_heigth(x, y);
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

uint32_t create_levelset(const int x, const int y, const int z, const int depth, const heightmap& heights, tree_t& tree)
{
	const int dim = 1 << depth;

	tree_t::node node;

	for(int _y = 0; _y < dim; ++_y)
		for (int _x = 0; _x < dim; ++_x)
			if (uint16_t height = heights.at(x + _x, y + _y); height >= z && height < z + dim)
				goto ACTIVE;

	return 0;

ACTIVE:

	const int half_dim = 1 << (depth - 1);

	if (depth != 1)//No leaf yet...
	{
		node.children[0] = create_levelset(x           , y           , z           , depth - 1, heights, tree);
		node.children[1] = create_levelset(x + half_dim, y           , z           , depth - 1, heights, tree);
		node.children[2] = create_levelset(x           , y + half_dim, z           , depth - 1, heights, tree);
		node.children[3] = create_levelset(x + half_dim, y + half_dim, z           , depth - 1, heights, tree);
		node.children[4] = create_levelset(x           , y           , z + half_dim, depth - 1, heights, tree);
		node.children[5] = create_levelset(x + half_dim, y           , z + half_dim, depth - 1, heights, tree);
		node.children[6] = create_levelset(x           , y + half_dim, z + half_dim, depth - 1, heights, tree);
		node.children[7] = create_levelset(x + half_dim, y + half_dim, z + half_dim, depth - 1, heights, tree);
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

uint32_t create_volume(const int x, const int y, const int z, const int depth, const heightmap& heights, tree_t& tree)
{
	const int dim = 1 << depth;

	tree_t::node node;

	for(int _y = 0; _y < dim; ++_y)
		for (int _x = 0; _x < dim; ++_x)
			if (z <= heights.at(x + _x, y + _y))
				goto ACTIVE;

	return 0;

ACTIVE:

	if (depth != 1)//No leaf yet...
	{
		const int half_dim = 1 << (depth - 1);

		node.children[0] = create_volume(x           , y           , z           , depth - 1, heights, tree);
		node.children[1] = create_volume(x + half_dim, y           , z           , depth - 1, heights, tree);
		node.children[2] = create_volume(x           , y + half_dim, z           , depth - 1, heights, tree);
		node.children[3] = create_volume(x + half_dim, y + half_dim, z           , depth - 1, heights, tree);
		node.children[4] = create_volume(x           , y           , z + half_dim, depth - 1, heights, tree);
		node.children[5] = create_volume(x + half_dim, y           , z + half_dim, depth - 1, heights, tree);
		node.children[6] = create_volume(x           , y + half_dim, z + half_dim, depth - 1, heights, tree);
		node.children[7] = create_volume(x + half_dim, y + half_dim, z + half_dim, depth - 1, heights, tree);
	}
	else//Leaf
	{
		node.children[0] = z     <= heights.at(x    , y    ) ? get_leaf_val() : 0;
		node.children[1] = z     <= heights.at(x + 1, y    ) ? get_leaf_val() : 0;
		node.children[2] = z     <= heights.at(x    , y + 1) ? get_leaf_val() : 0;
		node.children[3] = z     <= heights.at(x + 1, y + 1) ? get_leaf_val() : 0;
		node.children[4] = z + 1 <= heights.at(x    , y    ) ? get_leaf_val() : 0;
		node.children[5] = z + 1 <= heights.at(x + 1, y    ) ? get_leaf_val() : 0;
		node.children[6] = z + 1 <= heights.at(x    , y + 1) ? get_leaf_val() : 0;
		node.children[7] = z + 1 <= heights.at(x + 1, y + 1) ? get_leaf_val() : 0;
	}

	if (node.is_zero())
		printf("ERR");

	return tree.register_node(node);
}



template<typename Noise>
uint32_t fill_with(const int x, const int y, const int z, const int depth, tree_t& tree, Noise& noise)
{
	tree_t::node node;

	if (depth != 1)//No leaf yet...
	{
		const int half_dim = 1 << (depth - 1);

		node.children[0] = fill_with(x           , y           , z           , depth - 1, tree, noise);
		node.children[1] = fill_with(x + half_dim, y           , z           , depth - 1, tree, noise);
		node.children[2] = fill_with(x           , y + half_dim, z           , depth - 1, tree, noise);
		node.children[3] = fill_with(x + half_dim, y + half_dim, z           , depth - 1, tree, noise);
		node.children[4] = fill_with(x           , y           , z + half_dim, depth - 1, tree, noise);
		node.children[5] = fill_with(x + half_dim, y           , z + half_dim, depth - 1, tree, noise);
		node.children[6] = fill_with(x           , y + half_dim, z + half_dim, depth - 1, tree, noise);
		node.children[7] = fill_with(x + half_dim, y + half_dim, z + half_dim, depth - 1, tree, noise);
	}
	else//Leaf
	{
		node.children[0] = noise(x    , y    , z    );
		node.children[1] = noise(x + 1, y    , z    );
		node.children[2] = noise(x    , y + 1, z    );
		node.children[3] = noise(x + 1, y + 1, z    );
		node.children[4] = noise(x    , y    , z + 1);
		node.children[5] = noise(x + 1, y    , z + 1);
		node.children[6] = noise(x    , y + 1, z + 1);
		node.children[7] = noise(x + 1, y + 1, z + 1);
	}

	if (node.is_zero())
		return 0;

	return tree.register_node(node);
}

template<typename Noise>
void remove(tree_t& tree, Noise& noise)
{
	for (int z = 0; z < tree.dim; ++z)
		for (int y = 0; y < tree.dim; ++y)
			for (int x = 0; x < tree.dim; ++x)
				if (uint32_t val = noise(x, y, z); val == 0)
					tree.set(x, y, z, val);
}

struct splatter_noise
{
	float threshold;
	uint32_t active_v;
	uint32_t inactive_v = 0;
	float scale;
	OpenSimplexNoise m_noise;

	splatter_noise(float threshold, uint32_t active_v, uint32_t inactive_v, float scale = tree_t::voxel_dim, int64_t seed = 123456) : threshold{ threshold }, active_v{ active_v }, inactive_v{ inactive_v }, scale{ scale }, m_noise{ OpenSimplexNoise(seed) } {}

	uint32_t operator()(int x, int y, int z)
	{
		float fx = static_cast<float>(x) * scale;
		float fy = static_cast<float>(y) * scale;
		float fz = static_cast<float>(z) * scale;

		float val = static_cast<float>(noise(fx, fy, fz));

		return val >= threshold ? active_v : inactive_v;
	}
};

void initialize_h_octree(tree_t& tree)
{
	splatter_noise caverns(-0.6F, 1, 0, 1.0F / 64.0F, 1282);
	splatter_noise tunnels(-0.5F, 1, 0, 1.0F / 16.0F, 9767564);

	heightmap heights(tree.dim);
	
	tree.set_root(create_volume(0, 0, 0, tree.depth, heights, tree));

	for (uint16_t y = 0; y < tree.dim; ++y)
		for (uint16_t x = 0; x < tree.dim; ++x)
		{
			uint16_t z = heights.at(x, y);
			tree.set(x, y, z    , 2 + (std::rand() > RAND_MAX / 2));
			tree.set(x, y, z - 1, 4);
			tree.set(x, y, z - 2, 4);
		}

	//remove(tree, caverns);
	remove(tree, tunnels);
}



void test_och_h_octree()
{
	printf("\nEntered test_och_h_octree...\n");

	std::string voxel_loading_errmsg;

	och::voxel_data voxels("voxels.txt", voxel_loading_errmsg);

	if (voxel_loading_errmsg.size())
	{
		printf("\n%s\n", voxel_loading_errmsg.c_str());

		printf("\n...Exited test_och_hashed_octree\n");

		return;
	}

	printf("\nRead data for %i voxels\n", voxels.voxel_cnt);

	printf("\nConstructing tree\n");

	tree_t tree;

	printf("\nTree-depth: %i\n", tree.depth);

	printf("\nTree-dimension: %i\n", tree.dim);

	printf("\nTable-size: %s\n", och::abbreviate_byte_size(tree.table_bytes).c_str());

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	initialize_h_octree(tree);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	printf("\nTime taken to construct svodag:   %lli ms\n", std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count());

	printf("\nFilled entries in node-hashtable: %i of %i max\n", tree.get_fillcnt(), tree.table_capacity);

	printf("\nTotal node-count: %i\n", tree.get_nodecnt());

	printf("\nCompression-ratio: %f\n", (float)tree.get_nodecnt() / tree.get_fillcnt());

	printf("\nRequired size in memory: %s\n", och::abbreviate_byte_size(static_cast<size_t>(tree.get_fillcnt()) * 37).c_str());

	printf("\nUncompressed size: %s\n", och::abbreviate_byte_size(static_cast<size_t>(tree.get_nodecnt()) * 32).c_str());

	tree_window window(tree, voxels);

	if (window.Construct(screen_size_x, screen_size_y, pixel_size, pixel_size))
	{
		printf("\nOpening window...\n");

		while (ShowCursor(false) >= 0);

		window.Start();
	}
	else
		printf("\nCould not open window\n");

	printf("\n...Exited test_och_hashed_octree\n");
}
