#pragma once

#include <math.h>
#include "ofApp.h"

int IX(int x, int y)
{
	x = ofClamp(x, 0, ofApp::N - 1);
	y = ofClamp(y, 0, ofApp::N - 1);
	return x + y * ofApp::N;
}

struct Fluid
{
	int size;
	float dt;
	float diff;
	float visc;

	std::vector<float> s;
	std::vector<float> density;

	std::vector<float> Vx;
	std::vector<float> Vy;

	std::vector<float> Vx0;
	std::vector<float> Vy0;

	Fluid(float dt, float diffusion, float viscosity)
		: dt(dt), diff(diffusion), visc(viscosity), size(ofApp::N)
	{
		this->s = std::vector<float>(ofApp::N * ofApp::N);
		this->density = std::vector<float>(ofApp::N * ofApp::N);

		this->Vx = std::vector<float>(ofApp::N * ofApp::N);
		this->Vy = std::vector<float>(ofApp::N * ofApp::N);

		this->Vx0 = std::vector<float>(ofApp::N * ofApp::N);
		this->Vy0 = std::vector<float>(ofApp::N * ofApp::N);
	}

	void step()
	{
		diffuse(1, this->Vx0, this->Vx, this->visc, this->dt);
		diffuse(2, this->Vy0, this->Vy, this->visc, this->dt);

		project(this->Vx0, this->Vy0, this->Vx, this->Vy);

		advect(1, this->Vx, this->Vx0, this->Vx0, this->Vy0, this->dt);
		advect(2, this->Vy, this->Vy0, this->Vx0, this->Vy0, this->dt);

		project(this->Vx, this->Vy, this->Vx0, this->Vy0);

		diffuse(0, this->s, this->density, this->diff, this->dt);
		advect(0, this->density, this->s, this->Vx, this->Vy, this->dt);
	}

	void render_d()
	{
		for (int i = 0; i < ofApp::N; i++)
		{
			for (int j = 0; j < ofApp::N; j++)
			{
				float x = i * ofApp::SCALE;
				float y = j * ofApp::SCALE;
				float d = this->density[IX(i, j)];
				if (d != 0)
				{
					ofSetColor(ofColor::fromHsb(int(d+50)%255, 200, d));
					ofDrawRectangle(x, y, ofApp::SCALE, ofApp::SCALE);
				}
			}
		}
	}

	void fade_d()
	{
		for (int i = 0; i < this->density.size(); i++)
		{
			float d = this->density[i];
			this->density[i] = ofClamp(d - 0.02, 0, 255);
		}
	}

	void add_density(int x, int y, float amount)
	{
		this->density[IX(x, y)] += amount;
	}

	void add_velocity(int x, int y, float amount_x, float amount_y)
	{
		this->Vx[IX(x, y)] += amount_x;
		this->Vy[IX(x, y)] += amount_y;
	}

	void diffuse(int b, std::vector<float>& x, std::vector<float>& x0, float diff, float dt)
	{
		float a = dt * diff * (ofApp::N - 2) * (ofApp::N - 2);
		lin_solve(b, x, x0, a, 1 + 4 * a);
	}

	void project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& p, std::vector<float>& div)
	{
		for (int j = 1; j < ofApp::N - 1; j++) {
			for (int i = 1; i < ofApp::N - 1; i++) {
				div[IX(i, j)] = -0.5f * (
					velocX[IX(i + 1, j)]
					- velocX[IX(i - 1, j)]
					+ velocY[IX(i, j + 1)]
					- velocY[IX(i, j - 1)]
					) / ofApp::N;
				p[IX(i, j)] = 0;
			}
		}

		set_bnd(0, div);
		set_bnd(0, p);
		lin_solve(0, p, div, 1, 4);

		for (int j = 1; j < ofApp::N - 1; j++) {
			for (int i = 1; i < ofApp::N - 1; i++) {
				velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)]
					- p[IX(i - 1, j)]) * ofApp::N;
				velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)]
					- p[IX(i, j - 1)]) * ofApp::N;
			}
		}
		set_bnd(1, velocX);
		set_bnd(2, velocY);
	}

	void advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& velocX, std::vector<float>& velocY, float dt) {
		float i0, i1, j0, j1;

		float dtx = dt * (ofApp::N - 2);
		float dty = dt * (ofApp::N - 2);

		float s0, s1, t0, t1;
		float tmp1, tmp2, x, y;

		float Nfloat = ofApp::N;
		float ifloat, jfloat;
		int i, j;

		for (j = 1, jfloat = 1; j < ofApp::N - 1; j++, jfloat++) {
			for (i = 1, ifloat = 1; i < ofApp::N - 1; i++, ifloat++) {
				tmp1 = dtx * velocX[IX(i, j)];
				tmp2 = dty * velocY[IX(i, j)];
				x = ifloat - tmp1;
				y = jfloat - tmp2;

				if (x < 0.5f) x = 0.5f;
				if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
				i0 = floor(x);
				i1 = i0 + 1.0f;
				if (y < 0.5f) y = 0.5f;
				if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
				j0 = floor(y);
				j1 = j0 + 1.0f;

				s1 = x - i0;
				s0 = 1.0f - s1;
				t1 = y - j0;
				t0 = 1.0f - t1;

				int i0i = int(i0);
				int i1i = int(i1);
				int j0i = int(j0);
				int j1i = int(j1);

				// DOUBLE CHECK THIS!!!
				d[IX(i, j)] =
					s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
					s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
			}
		}

		set_bnd(b, d);
	}

	void lin_solve(int b, std::vector<float>& x, std::vector<float>& x0, float a, float c)
	{
		float cRecip = 1.0 / c;
		for (int k = 0; k < ofApp::iter; k++)
		{
			for (int j = 1; j < ofApp::N - 1; j++)
			{
				for (int i = 1; i < ofApp::N - 1; i++)
				{
					x[IX(i, j)] =
						(x0[IX(i, j)]
							+ a * (x[IX(i + 1, j)]
								+ x[IX(i - 1, j)]
								+ x[IX(i, j + 1)]
								+ x[IX(i, j - 1)]
								)) * cRecip;
				}
			}

			set_bnd(b, x);
		}
	}

	void set_bnd(int b, std::vector<float>& x)
	{
		for (int i = 1; i < ofApp::N - 1; i++)
		{
			x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
			x[IX(i, ofApp::N - 1)] = b == 2 ? -x[IX(i, ofApp::N - 2)] : x[IX(i, ofApp::N - 2)];
		}

		for (int j = 1; j < ofApp::N - 1; j++)
		{
			x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
			x[IX(ofApp::N - 1, j)] = b == 1 ? -x[IX(ofApp::N - 2, j)] : x[IX(ofApp::N - 2, j)];
		}

		x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
		x[IX(0, ofApp::N - 1)] = 0.5f * (x[IX(1, ofApp::N - 1)] + x[IX(0, ofApp::N - 2)]);
		x[IX(ofApp::N - 1, 0)] = 0.5f * (x[IX(ofApp::N - 2, 0)] + x[IX(ofApp::N - 1, 1)]);
		x[IX(ofApp::N - 1, ofApp::N - 1)] = 0.5f * (x[IX(ofApp::N - 2, ofApp::N - 1)] + x[IX(ofApp::N - 1, ofApp::N - 2)]);
	}
};