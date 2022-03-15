#pragma once
#include <math.h>
#include <Eigen/Dense>


namespace quatmath_consts
{
    const float eps = 1e-6f;
    const float pi = 3.1415f;
}

Eigen::Vector4f quatFromEul(const Eigen::Vector3f& eul)
{
	 float roll = eul[0];
	 float pitch = eul[1];
	 float yaw = eul[2];

	 float cy = cos(yaw * 0.5f);
	 float sy = sin(yaw * 0.5f);
	 float cr = cos(roll * 0.5f);
	 float sr = sin(roll * 0.5f);
	 float cp = cos(pitch * 0.5f);
	 float sp = sin(pitch * 0.5f);

	 float q0 = cy * cr * cp + sy * sr * sp;
	 float q1 = cy * sr * cp - sy * cr * sp;
	 float q2 = cy * cr * sp + sy * sr * cp;
	 float q3 = sy * cr * cp - cy * sr * sp;

	 return Eigen::Vector4f(q0, q1, q2 ,q3);
}

Eigen::Vector4f quatMultiply(const Eigen::Vector4f& q, const Eigen::Vector4f& r)
{
	Eigen::Vector4f p;
	p[0] = r[0] * q[0] - r[1] * q[1] - r[2] * q[2] - r[3] * q[3];
	p[1] = r[0] * q[1] + r[1] * q[0] - r[2] * q[3] + r[3] * q[2];
	p[2] = r[0] * q[2] + r[1] * q[3] + r[2] * q[0] - r[3] * q[1];
	p[3] = r[0] * q[3] - r[1] * q[2] + r[2] * q[1] + r[3] * q[0];
	return p;
}

Eigen::Vector4f quatInverse(const Eigen::Vector4f& q)
{
	Eigen::Vector4f qDual (q[0], -q[1], -q[2], -q[3]);
	return qDual;
}

Eigen::Vector3f quatRotate(const Eigen::Vector4f& q, const Eigen::Vector3f& v)
{
	Eigen::Vector4f qv(0.0f, v[0], v[1], v[2]);

	Eigen::Vector4f qDual = quatInverse(q);
	Eigen::Vector4f qv1 = quatMultiply(qv, qDual);
	Eigen::Vector4f qv2 = quatMultiply(q, qv1);

	return Eigen::Vector3f(qv2[1], qv2[2], qv2[3]);
}

Eigen::Matrix<float, 3, 3> quatToMatrix(const Eigen::Vector4f& q)
{
	float R11 = 1.f - 2.f * q(2) * q(2) - 2.f * q(3) * q(3);
	float R12 = 2.f * q(1) * q(2) - 2.f * q(3) * q(0);
	float R13 = 2.f * q(1) * q(3) + 2.f * q(2) * q(0);

	float R21 = 2.f * q(1) * q(2) + 2.f * q(3) * q(0);
	float R22 = 1.f - 2.f * q(1) * q(1) - 2.f * q(3) * q(3);
	float R23 = 2.f * q(2) * q(3) - 2.f * q(1) * q(0);

	float R31 = 2.f * q(1) * q(3) - 2.f * q(2) * q(0);
	float R32 = 2.f * q(2) * q(3) + 2.f * q(1) * q(0);
	float R33 = 1.f - 2.f * q(1) * q(1) - 2.f * q(2) * q(2);

	Eigen::Matrix<float, 3, 3> R;
	R << R11, R12, R13, R21, R22, R23, R31, R32, R33;
	return R;
}

float shortestRotation(float from, float to)
{
    float da = fmod((to - from + quatmath_consts::pi), 2*quatmath_consts::pi) - quatmath_consts::pi;
    if (da < -quatmath_consts::pi)
    {
        da = da + 2*quatmath_consts::pi;
    }
    return da;
}

Eigen::Vector4f wgsToEnuQuat(float lat, float lon)
{
    double pi = 3.1415;
    Eigen::Vector4f qlat = quatFromEul(Eigen::Vector3f(pi/2 - lat, 0, 0));
    Eigen::Vector4f qlon = quatFromEul(Eigen::Vector3f(0, 0, pi/2 + lon));
                 
    Eigen::Vector4f res = quatMultiply(qlon, qlat);
    
    return res;
}
