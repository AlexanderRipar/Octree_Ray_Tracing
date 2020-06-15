#include <iostream>
#include <cstdio>
#include <immintrin.h>

#include "och_float.h"
#include "test_och_h_octree.h"
#include "test_och_octree.h"
#include "och_string_util.h"

/*constexpr int t_rep = 1;

void bresenham_own(pixel_grid& p, float x0, float y0, float x1, float y1) {

	int x = static_cast<int>(x0);
	int y = static_cast<int>(y0);

	int d_x = static_cast<int>(x1 - x0);
	int d_y = static_cast<int>(y1 - y0);

	float d_err = (y1 - y0) / (x1 - x0);

	float err = y0 - y + d_err * (x - x0);

	while (x <= x1) {
		
		err += d_err;

		if (err >= 1) {
			p.set(x, y++, t_rep);
			err--;
		}

		p.set(x++, y, t_rep);
	}

}*/

/*void bresenham_test(pixel_grid& p, float x1, float y1, float x2, float y2) {

	int i;               // loop counter
	int ystep, xstep;    // the step on y and x axis
	int error;           // the error accumulated during the increment
	int errorprev;       // *vision the previous value of the error variable
	int y = static_cast<int>(y1), x = static_cast<int>(x1);  // the line points
	int ddy, ddx;        // compulsory variables: the double values of dy and dx
	int dx = static_cast<int>(x2 - x1);
	int dy = static_cast<int>(y2 - y1);
	p.set(x, y, t_rep);  // first point
	// NB the last point can't be here, because of its previous point (which has to be verified)
	if (dy < 0) {
		ystep = -1;
		dy = -dy;
	}
	else
		ystep = 1;

	if (dx < 0) {
		xstep = -1;
		dx = -dx;
	}
	else
		xstep = 1;

	ddy = 2 * dy;  // work with double values for full precision
	ddx = 2 * dx;

	if (ddx >= ddy) {  // first octant (0 <= slope <= 1)

		// compulsory initialization (even for errorprev, needed when dx==dy)
		errorprev = error = dx;  // start in the middle of the square

		for (i = 0; i < dx; i++) {  // do not use the first point (already done)

			x += xstep;
			error += ddy;

			if (error > ddx) {  // increment y if AFTER the middle ( > )

				y += ystep;
				error -= ddx;

				// three cases (octant == right->right-top for directions below):
				if (error + errorprev < ddx)		// bottom square also
					p.set(x, y - ystep, t_rep);
				else if (error + errorprev > ddx)	// left square also
					p.set(x - xstep, y, t_rep);
				else {								// corner: bottom and left squares also
					p.set(x, y - ystep, t_rep);
					p.set(x - xstep, y, t_rep);
				}
			}

			p.set(x, y, t_rep);
			errorprev = error;
		}
	}
	else {  // the same as above

		errorprev = error = dy;

		for (i = 0; i < dy; i++) {

			y += ystep;
			error += ddx;

			if (error > ddy) {

				x += xstep;
				error -= ddy;

				if (error + errorprev < ddy)
					p.set(x - xstep, y, t_rep);
				else if (error + errorprev > ddy)
					p.set(x, y - ystep, t_rep);
				else {
					p.set(x - xstep, y, t_rep);
					p.set(x, y - ystep, t_rep);
				}
			}

			p.set(x, y, t_rep);
			errorprev = error;
		}
	}
	_ASSERT((y == y2) && (x == x2));  // the last point (y2,x2) has to be the same with the last point of the algorithm
}*/

/*
std::vector<int> UpTo(int n, int offset = 0)
{
	std::vector<int> retval(n);
	for (int ii = 0; ii < n; ++ii)
		retval[ii] = ii + offset;
	return retval;
}

struct JohnsonTrotterState_
{
	std::vector<int> values_;
	std::vector<int> positions_;	// size is n+1, first element is not used
	std::vector<bool> directions_;
	int sign_;

	JohnsonTrotterState_(int n) : values_(UpTo(n, 1)), positions_(UpTo(n + 1, -1)), directions_(n + 1, false), sign_(1) {}

	int LargestMobile() const	// returns 0 if no mobile integer exists
	{
		for (int r = values_.size(); r > 0; --r)
		{
			const int loc = positions_[r] + (directions_[r] ? 1 : -1);
			if (loc >= 0 && loc < values_.size() && values_[loc] < r)
				return r;
		}
		return 0;
	}

	bool IsComplete() const { return LargestMobile() == 0; }

	void operator++()	// implement Johnson-Trotter algorithm
	{
		const int r = LargestMobile();
		const int rLoc = positions_[r];
		const int lLoc = rLoc + (directions_[r] ? 1 : -1);
		const int l = values_[lLoc];
		// do the swap
		std::swap(values_[lLoc], values_[rLoc]);
		std::swap(positions_[l], positions_[r]);
		sign_ = -sign_;
		// change directions
		for (auto pd = directions_.begin() + r + 1; pd != directions_.end(); ++pd)
			*pd = !*pd;
	}
};

void permutation_solver()
{
#define FISH	(s.values_[0] - 1)
#define BALL	(s.values_[1] - 1)
#define HEART	(s.values_[2] - 1)
#define FIRE	(s.values_[3] - 1)
#define PARTY	(s.values_[4] - 1)
#define POOP	(s.values_[5] - 1)
#define MEAT	(s.values_[6] - 1)
#define PERSON	(s.values_[7] - 1)
#define CAMERA	(s.values_[8] - 1)
#define SMILEY	(s.values_[9] - 1)

	int no_of_solutions = 0;

	std::vector<int> sol(10);

	JohnsonTrotterState_ s(10);

	std::cout << "\nChecking permutations...\n";

	do
	{
		if ((FIRE + 10 * HEART + 100 * FISH + 1000 * BALL + 10000 * FISH) - (PERSON + 10 * MEAT + 100 * POOP + 1000 * PARTY) == SMILEY + 10 * BALL + 100 * FISH + 1000 * CAMERA)
			if ((PARTY + 100 * BALL + 1000 * MEAT) - (HEART + 100 * SMILEY) == PERSON + 10 * CAMERA + 100 * FIRE + 1000 * FISH)
			{
				++no_of_solutions;
				std::cout << "\nFISH:   " << FISH << "\nBALL:   " << BALL << "\nHEART:  " << HEART << "\nFIRE:   " << 
					FIRE << "\nPARTY:  " << PARTY << "\nPOOP:   " << POOP << "\nMEAT:   " << MEAT << "\nPERSON: " << PERSON << "\nCAMERA: " << CAMERA << "\nSMILEY: " << SMILEY << '\n';

				if (no_of_solutions == 1)
					sol = s.values_;
			}

		++s;
	} while (!s.IsComplete());

	std::cout << "\nFinished permutations...\n";

	if (no_of_solutions == 1)
	{
		std::cout << "\nUnique solution found! The formulae are\n\n\t";

		for (auto& i : sol)
			--i;

		int first_a = sol[3] + 10 * sol[2] + 100 * sol[0] + 1000 * sol[1] + 10000 * sol[0];
		int first_b = sol[7] + 10 * sol[6] + 100 * sol[5] + 1000 * sol[4];
		int first_r = sol[9] + 10 * sol[1] + 100 * sol[0] + 1000 * sol[8];

		int second_a = sol[4] + 10 * sol[5] + 100 * sol[1] + 1000 * sol[6];
		int second_b = sol[2] + 10 * sol[5] + 100 * sol[9];
		int second_r = sol[7] + 10 * sol[8] + 100 * sol[3] + 1000 * sol[0];

		std::cout << first_a << " - " << first_b << " = " << first_r << "\n\nand\n\n\t" << second_a << " - " << second_b << " = " << second_r << "\n\n";
	}

#undef FISH
#undef BALL
#undef HEART
#undef FIRE
#undef PARTY
#undef POOP
#undef MEAT
#undef PERSON
#undef CAMERA
#undef SMILEY
}
*/

/*
void find_words()
{
	const std::string grid_str = 
		"AWAKIIHTOHSLDQTEAUSP"
		"YGSSDFCOFIREGUEKTLGU"
		"ASAOLZRFSERILWTHSIUS"
		"LBCZJELTFAIDDTLNIHPI"
		"TAOSBLMCSIHCADOONJNZ"
		"PLAAICHSSESDOOIPAOEU"
		"INWBDTIDWSMSAEHEYEWE"
		"OSRIZWSEIHANDERSSWAO"
		"LAOAHALCHITLTFFHVOYE"
		"EITEQKTNIVHBTNISTAAB"
		"EOMTLAIHUDSEIBOSDDSR"
		"WEMEHITCINNVLNUSCNUE"
		"AJGAIZHANEBTAATCSAEI"
		"UGGLTADIUSSIMBLTLHIH"
		"UDRRLSIOKNRESIEQSIOV"
		"PEESTHENBAYOCHNTOPOE"
		"IAETUILEXECDLIEEVODS"
		"AUPIODOEGREBMEALFSSE"
		"TLMIHYOUMCCGDLGLNTAE"
		"UIEARDODRSSEVTPSOZLE";

	char grid[20][20];

	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < 20; ++j)
			grid[i][j] = grid_str[i * 20 + j];

	std::vector<std::string> lines;

	//STRAIGHT
	for (int i = 0; i < 20; ++i)
	{
		std::string h_fwd;
		std::string h_bwd;
		std::string v_dwd;
		std::string v_uwd;

		for (int j = 0; j < 20; ++j)
		{
			h_fwd.push_back(grid[i][j]);
			h_bwd.push_back(grid[i][19 - j]);
			v_dwd.push_back(grid[j][i]);
			v_uwd.push_back(grid[19 - j][i]);

		}

			lines.push_back(h_fwd);
			lines.push_back(h_bwd);
			lines.push_back(v_dwd);
			lines.push_back(v_uwd);
	}

	//DIAGONAL
	for (int i = 0; i < 20; ++i)
	{
		for (int j = 0; j < 20; ++j);//TODO
	}

	while (1)
	{
		std::cout << "Enter word ";

		std::string to_find;
		std::getline(std::cin, to_find);

		for (auto& c : to_find)
			c = toupper(c);

		if (to_find == "EXIT")
			return;
		else if (to_find == "")
			continue;
		else if (to_find.length() > 20)
		{
			std::cout << "Word too long\n\n";
			continue;
		}

		for (auto& l : lines)
		{
			for (int i = 0; i < 20 - to_find.size() - 1; ++i)
			{
				bool is_equal = true;


				for (int c = 0; c < to_find.size(); ++c)
				{
					if (l[i + c] != to_find[c])
					{
						if (i + c + 1 < 20 && l[i + c + 1] == to_find[c])
						{
							++i;
						}
						else
						{
							is_equal = false;
							break;
						}
					}
				}

				if (is_equal)
				{
					std::string eq = l.substr(i - 1, to_find.size() + 1);

					std::cout << "Match FOUND: " << eq << "\n\n";

					goto END_SEARCH;
				}
			}
		}

		std::cout << "No match found\n\n";

	END_SEARCH:;
	}
}

std::string trace_word(const std::string& to_find, int x, int y, int inc_x, int inc_y, const char grid[20][20])
{
	std::string ret;

	ret.push_back(grid[x][y]);

	int _x = x + inc_x, _y = y + inc_y;

	bool has_skipped = false;

	int pos = 1;

	while (pos < to_find.size() && _x < 20 && _x >= 0 && _y < 20 && _y >= 0)
	{
		if (grid[_x][_y] != to_find[pos])
			if (grid[_x + inc_x][_y + inc_y] != to_find[pos] || has_skipped)
				return "";
			else
			{
				has_skipped = true;
				ret.push_back(grid[_x][_y]);
				_x += inc_x;
				_y += inc_y;
			}

		ret.push_back(grid[_x][_y]);

		++pos;
		_x += inc_x;
		_y += inc_y;
	}

	return pos == to_find.size() ? ret : "";
}

void find_words_2()
{
	const std::string grid_str =
		"AWAKIIHTOHSLDQTEAUSP"
		"YGSSDFCOFIREGUEKTLGU"
		"ASAOLZRFSERILWTHSIUS"
		"LBCZJELTFAIDDTLNIHPI"
		"TAOSBLMCSIHCADOONJNZ"
		"PLAAICHSSESDOOIPAOEU"
		"INWBDTIDWSMSAEHEYEWE"
		"OSRIZWSEIHANDERSSWAO"
		"LAOAHALCHITLTFFHVOYE"
		"EITEQKTNIVHBTNISTAAB"
		"EOMTLAIHUDSEIBOSDDSR"
		"WEMEHITCINNVLNUSCNUE"
		"AJGAIZHANEBTAATCSAEI"
		"UGGLTADIUSSIMBLTLHIH"
		"UDRRLSIOKNRESIEQSIOV"
		"PEESTHENBAYOCHNTOPOE"
		"IAETUILEXECDLIEEVODS"
		"AUPIODOEGREBMEALFSSE"
		"TLMIHYOUMCCGDLGLNTAE"
		"UIEARDODRSSEVTPSOZLE";

	char grid[20][20];

	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < 20; ++j)
			grid[i][j] = grid_str[i * 20 + j];

	while (1)
	{
		std::cout << "Enter word ";

		std::string to_find;
		std::getline(std::cin, to_find);

		for (auto& c : to_find)
			c = toupper(c);

		if (to_find == "EXIT")
			return;
		else if (to_find == "")
			continue;
		else if (to_find.length() > 20)
		{
			std::cout << "Word too long\n\n";
			continue;
		}

		int len = 0;
		int init_idx_x;
		int init_idx_y;
		
		enum class dir
		{
			n, ne, e, se, s, sw, w, nw
		};

		dir d;

		//For every field
			//if first letter matches
				//for every direction
					//if not outside
						//if second letter matches
							//match rest

		for (int i = 0; i < 20; ++i)
			for (int j = 0; j < 20; ++j)
			{
				if (grid[i][j] == to_find[0])
				{
					std::string n = trace_word(to_find, i, j, 0, -1, grid);
					if (n != "")
					{
						std::cout << "MATCH FOUND: " << n << "\n\n";
						goto SEARCH_END;
					}

					std::string ne = trace_word(to_find, i, j, 1, -1, grid);
					if (ne != "")
					{
						std::cout << "MATCH FOUND: " << ne << "\n\n";
						goto SEARCH_END;
					}

					std::string e = trace_word(to_find, i, j, 1, 0, grid);
					if (e != "")
					{
						std::cout << "MATCH FOUND: " << e << "\n\n";
						goto SEARCH_END;
					}

					std::string se = trace_word(to_find, i, j, 1, 1, grid);
					if (se != "")
					{
						std::cout << "MATCH FOUND: " << se << "\n\n";
						goto SEARCH_END;
					}

					std::string s = trace_word(to_find, i, j, 0, 1, grid);
					if (s != "")
					{
						std::cout << "MATCH FOUND: " << s << "\n\n";
						goto SEARCH_END;
					}

					std::string sw = trace_word(to_find, i, j, -1, 1, grid);
					if (sw != "")
					{
						std::cout << "MATCH FOUND: " << sw << "\n\n";
						goto SEARCH_END;
					}

					std::string w = trace_word(to_find, i, j, -1, 0, grid);
					if (w != "")
					{
						std::cout << "MATCH FOUND: " << w << "\n\n";
						goto SEARCH_END;
					}

					std::string nw = trace_word(to_find, i, j, -1, -1, grid);
					if (nw != "")
					{
						std::cout << "MATCH FOUND: " << nw << "\n\n";
						goto SEARCH_END;
					}
				}
			}

		std::cout << "No match found\n\n";

	SEARCH_END:;
	}
}
*/

int main()
{
	test_och_hashed_octree();

	//test_och_octree();





















	/*
	constexpr int item_count = 7;
	
	och::ocd_item items[item_count];
	
	items[0].name = "zassgsgf";
	items[1].name = "asdad";
	items[2].name = "aaa";
	items[3].name = "bsfg";
	items[4].name = "bsfgg";
	items[5].name = "asdf";
	items[6].name = "zassgsgf";
	
	och::write_ocd_file("a", items, item_count);
	*/

	/*
	#include <spasdt.h>

#define staindetect true

	bool staindetect = true;

	std::cout << "Lixi, are you a stain?\n" << std::boolalpha;
	std::string Lixi = "spasdt";
	if (Lixi == "spasdt") std::cout << staindetect << '\n';
	*/
}
