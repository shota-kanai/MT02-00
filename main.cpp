#define _USE_MATH_DEFINES
#include <math.h>
#include <Novice.h>
#include <assert.h>
#include <cmath>
#include<Vector3.h>
#include <imgui.h>
#include <algorithm>
#include<numbers>
const char kWindowTitle[] = "GC1C_05_カナイショウタ_タイトル";

static const int kRowHeight = 20;
static const int kColumnWidth = 60;

static const int kWindowWith = 1280;
static const int kWindowHigat = 720;


//struct Vector3 {
//	float x, y, z;
//};

struct Matrix3x3
{
	float m[3][3];
};

struct Matrix4x4
{
	float m[4][4];
};


struct Segment
{
	Vector3 origin;
	Vector3 diff;
};

struct Line
{
	Vector3 origin;
	Vector3 diff;
};

struct Ray
{
	Vector3 origin;
	Vector3 diff;
};
struct Sphere
{
	Vector3 center;
	float radius;
};
struct  Plane
{
	Vector3 normal;
	float distance;
};
struct Triangle
{
	Vector3 vertices[3];
};
struct AABB
{
	Vector3 min;
	Vector3 max;
};

Vector3 operator+(Vector3 num1, Vector3 num2) {
	num1.x += num2.x;
	num1.y += num2.y;
	num1.z += num2.z;

	return num1;
}

Vector3 operator-(Vector3 num1, Vector3 num2) {
	num1.x -= num2.x;
	num1.y -= num2.y;
	num1.z -= num2.z;

	return num1;
}

void VectorScreenPrintf(int x, int y, const Vector3& vector, const char* label) {
	Novice::ScreenPrintf(x, y, "%.02f", vector.x);
	Novice::ScreenPrintf(x + kColumnWidth, y, "%.02f", vector.y);
	Novice::ScreenPrintf(x + kColumnWidth * 2, y, "%.02f", vector.z);
	Novice::ScreenPrintf(x + kColumnWidth * 3, y, "%s", label);

}

void MatrixScreenPrintf(int x, int y, const Matrix4x4& matrix, const char* label) {

	Novice::ScreenPrintf(x, y, "%s", label);
	for (int row = 0; row < 4; ++row) {
		for (int column = 0; column < 4; ++column) {
			Novice::ScreenPrintf(x + column * kColumnWidth, (y + row * kRowHeight) + 20, "%6.02f", matrix.m[row][column]);
		}
	}
}



//加算
Vector3 Add(const Vector3& v1, const Vector3& v2) {
	return Vector3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);

}



//減算
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	return Vector3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}



//スカラー倍
Vector3 Multiply(float scalar,const Vector3& v) {

	Vector3 result;
	result.x = scalar * v.x;
	result.y = scalar * v.y;
	result.z = scalar * v.z;
	return result;

}

//内積
float Dot(const Vector3& v1, const Vector3& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

//長さ
float Length(const Vector3& v) {
	return sqrtf(Dot(v, v));
}


//正規化
inline Vector3 Normalize(const Vector3& v) {
	float v2 = 0.0f;
	v2 = Length(v);
	return{ v.x / v2,v.y / v2,v.z / v2 };
}


//行列の加法
Matrix4x4 Add(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m[i][j] = m1.m[i][j] + m2.m[i][j];
		}
	}

	return result;
}

//行列の減法
Matrix4x4 Subtract(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m[i][j] = m1.m[i][j] - m2.m[i][j];
		}
	}

	return result;
}

//行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			result.m[i][j] = 0.0f;
			for (int k = 0; k < 4; ++k) {
				result.m[i][j] += m1.m[i][k] * m2.m[k][j];
			}
		}
	}

	return result;
}

// 逆行列
Matrix4x4 Invers(const Matrix4x4& m)
{
	Matrix4x4 result{};
	float determinant = m.m[0][0] * (m.m[1][1] * m.m[2][2] * m.m[3][3] +
		m.m[2][1] * m.m[3][2] * m.m[1][3] +
		m.m[3][1] * m.m[1][2] * m.m[2][3] -
		m.m[3][1] * m.m[2][2] * m.m[1][3] -
		m.m[2][1] * m.m[1][2] * m.m[3][3] -
		m.m[1][1] * m.m[3][2] * m.m[2][3]) -
		m.m[0][1] * (m.m[1][0] * m.m[2][2] * m.m[3][3] +
			m.m[2][0] * m.m[3][2] * m.m[1][3] +
			m.m[3][0] * m.m[1][2] * m.m[2][3] -
			m.m[3][0] * m.m[2][2] * m.m[1][3] -
			m.m[2][0] * m.m[1][2] * m.m[3][3] -
			m.m[1][0] * m.m[3][2] * m.m[2][3]) +
		m.m[0][2] * (m.m[1][0] * m.m[2][1] * m.m[3][3] +
			m.m[2][0] * m.m[3][1] * m.m[1][3] +
			m.m[3][0] * m.m[1][1] * m.m[2][3] -
			m.m[3][0] * m.m[2][1] * m.m[1][3] -
			m.m[2][0] * m.m[1][1] * m.m[3][3] -
			m.m[1][0] * m.m[3][1] * m.m[2][3]) -
		m.m[0][3] * (m.m[1][0] * m.m[2][1] * m.m[3][2] +
			m.m[2][0] * m.m[3][1] * m.m[1][2] +
			m.m[3][0] * m.m[1][1] * m.m[2][2] -
			m.m[3][0] * m.m[2][1] * m.m[1][2] -
			m.m[2][0] * m.m[1][1] * m.m[3][2] -
			m.m[1][0] * m.m[3][1] * m.m[2][2]);



	if (determinant != 0) {
		result.m[0][0] = (m.m[1][1] * m.m[2][2] * m.m[3][3] +
			m.m[2][1] * m.m[3][2] * m.m[1][3] +
			m.m[3][1] * m.m[1][2] * m.m[2][3] -
			m.m[3][1] * m.m[2][2] * m.m[1][3] -
			m.m[2][1] * m.m[1][2] * m.m[3][3] -
			m.m[1][1] * m.m[3][2] * m.m[2][3]) /
			determinant;

		result.m[0][1] = -(m.m[0][1] * m.m[2][2] * m.m[3][3] +
			m.m[2][1] * m.m[3][2] * m.m[0][3] +
			m.m[3][1] * m.m[0][2] * m.m[2][3] -
			m.m[3][1] * m.m[2][2] * m.m[0][3] -
			m.m[2][1] * m.m[0][2] * m.m[3][3] -
			m.m[0][1] * m.m[3][2] * m.m[2][3]) /
			determinant;

		result.m[0][2] = (m.m[0][1] * m.m[1][2] * m.m[3][3] +
			m.m[1][1] * m.m[3][2] * m.m[0][3] +
			m.m[3][1] * m.m[0][2] * m.m[1][3] -
			m.m[3][1] * m.m[1][2] * m.m[0][3] -
			m.m[1][1] * m.m[0][2] * m.m[3][3] -
			m.m[0][1] * m.m[3][2] * m.m[1][3]) /
			determinant;

		result.m[0][3] = -(m.m[0][1] * m.m[1][2] * m.m[2][3] +
			m.m[1][1] * m.m[2][2] * m.m[0][3] +
			m.m[2][1] * m.m[0][2] * m.m[1][3] -
			m.m[2][1] * m.m[1][2] * m.m[0][3] -
			m.m[1][1] * m.m[0][2] * m.m[2][3] -
			m.m[0][1] * m.m[2][2] * m.m[1][3]) /
			determinant;


		result.m[1][0] = -(m.m[1][0] * m.m[2][2] * m.m[3][3] +
			m.m[2][0] * m.m[3][2] * m.m[1][3] +
			m.m[3][0] * m.m[1][2] * m.m[2][3] -
			m.m[3][0] * m.m[2][2] * m.m[1][3] -
			m.m[2][0] * m.m[1][2] * m.m[3][3] -
			m.m[1][0] * m.m[3][2] * m.m[2][3]) /
			determinant;

		result.m[1][1] = (m.m[0][0] * m.m[2][2] * m.m[3][3] +
			m.m[2][0] * m.m[3][2] * m.m[0][3] +
			m.m[3][0] * m.m[0][2] * m.m[2][3] -
			m.m[3][0] * m.m[2][2] * m.m[0][3] -
			m.m[2][0] * m.m[0][2] * m.m[3][3] -
			m.m[0][0] * m.m[3][2] * m.m[2][3]) /
			determinant;

		result.m[1][2] = -(m.m[0][0] * m.m[1][2] * m.m[3][3] +
			m.m[1][0] * m.m[3][2] * m.m[0][3] +
			m.m[3][0] * m.m[0][2] * m.m[1][3] -
			m.m[3][0] * m.m[1][2] * m.m[0][3] -
			m.m[1][0] * m.m[0][2] * m.m[3][3] -
			m.m[0][0] * m.m[3][2] * m.m[1][3]) /
			determinant;

		result.m[1][3] = (m.m[0][0] * m.m[1][2] * m.m[2][3] +
			m.m[1][0] * m.m[2][2] * m.m[0][3] +
			m.m[2][0] * m.m[0][2] * m.m[1][3] -
			m.m[2][0] * m.m[1][2] * m.m[0][3] -
			m.m[1][0] * m.m[0][2] * m.m[2][3] -
			m.m[0][0] * m.m[2][2] * m.m[1][3]) /
			determinant;


		result.m[2][0] = (m.m[1][0] * m.m[2][1] * m.m[3][3] +
			m.m[2][0] * m.m[3][1] * m.m[1][3] +
			m.m[3][0] * m.m[1][1] * m.m[2][3] -
			m.m[3][0] * m.m[2][1] * m.m[1][3] -
			m.m[2][0] * m.m[1][1] * m.m[3][3] -
			m.m[1][0] * m.m[3][1] * m.m[2][3]) /
			determinant;

		result.m[2][1] = -(m.m[0][0] * m.m[2][1] * m.m[3][3] +
			m.m[2][0] * m.m[3][1] * m.m[0][3] +
			m.m[3][0] * m.m[0][1] * m.m[2][3] -
			m.m[3][0] * m.m[2][1] * m.m[0][3] -
			m.m[2][0] * m.m[0][1] * m.m[3][3] -
			m.m[0][0] * m.m[3][1] * m.m[2][3]) /
			determinant;

		result.m[2][2] = (m.m[0][0] * m.m[1][1] * m.m[3][3] +
			m.m[1][0] * m.m[3][1] * m.m[0][3] +
			m.m[3][0] * m.m[0][1] * m.m[1][3] -
			m.m[3][0] * m.m[1][1] * m.m[0][3] -
			m.m[1][0] * m.m[0][1] * m.m[3][3] -
			m.m[0][0] * m.m[3][1] * m.m[1][3]) /
			determinant;

		result.m[2][3] = -(m.m[0][0] * m.m[1][1] * m.m[2][3] +
			m.m[1][0] * m.m[2][1] * m.m[0][3] +
			m.m[2][0] * m.m[0][1] * m.m[1][3] -
			m.m[2][0] * m.m[1][1] * m.m[0][3] -
			m.m[1][0] * m.m[0][1] * m.m[2][3] -
			m.m[0][0] * m.m[2][1] * m.m[1][3]) /
			determinant;

		result.m[3][0] = -(m.m[1][0] * m.m[2][1] * m.m[3][2] +
			m.m[2][0] * m.m[3][1] * m.m[1][2] +
			m.m[3][0] * m.m[1][1] * m.m[2][2] -
			m.m[3][0] * m.m[2][1] * m.m[1][2] -
			m.m[2][0] * m.m[1][1] * m.m[3][2] -
			m.m[1][0] * m.m[3][1] * m.m[2][2]) /
			determinant;

		result.m[3][1] = (m.m[0][0] * m.m[2][1] * m.m[3][2] +
			m.m[2][0] * m.m[3][1] * m.m[0][2] +
			m.m[3][0] * m.m[0][1] * m.m[2][2] -
			m.m[3][0] * m.m[2][1] * m.m[0][2] -
			m.m[2][0] * m.m[0][1] * m.m[3][2] -
			m.m[0][0] * m.m[3][1] * m.m[2][2]) /
			determinant;

		result.m[3][2] = -(m.m[0][0] * m.m[1][1] * m.m[3][2] +
			m.m[1][0] * m.m[3][1] * m.m[0][2] +
			m.m[3][0] * m.m[0][1] * m.m[1][2] -
			m.m[3][0] * m.m[1][1] * m.m[0][2] -
			m.m[1][0] * m.m[0][1] * m.m[3][2] -
			m.m[0][0] * m.m[3][1] * m.m[1][2]) /
			determinant;

		result.m[3][3] = (m.m[0][0] * m.m[1][1] * m.m[2][2] +
			m.m[1][0] * m.m[2][1] * m.m[0][2] +
			m.m[2][0] * m.m[0][1] * m.m[1][2] -
			m.m[2][0] * m.m[1][1] * m.m[0][2] -
			m.m[1][0] * m.m[0][1] * m.m[2][2] -
			m.m[0][0] * m.m[2][1] * m.m[1][2]) /
			determinant;
	}

	return result;
}


//転置行列
Matrix4x4 Transpose(const Matrix4x4& m) {

	Matrix4x4 result = {};

	result.m[0][0] = m.m[0][0];
	result.m[0][1] = m.m[1][0];
	result.m[0][2] = m.m[2][0];
	result.m[0][3] = m.m[3][0];

	result.m[1][0] = m.m[0][1];
	result.m[1][1] = m.m[1][1];
	result.m[1][2] = m.m[2][1];
	result.m[1][3] = m.m[3][1];

	result.m[2][0] = m.m[0][2];
	result.m[2][1] = m.m[1][2];
	result.m[2][2] = m.m[2][2];
	result.m[2][3] = m.m[3][2];

	result.m[3][0] = m.m[0][3];
	result.m[3][1] = m.m[1][3];
	result.m[3][2] = m.m[2][3];
	result.m[3][3] = m.m[3][3];

	return result;
}

//単位行列
Matrix4x4 MakeIdentity4x4() {

	Matrix4x4 result = {};

	result.m[0][0] = 1;
	result.m[0][1] = 0;
	result.m[0][2] = 0;
	result.m[0][3] = 0;
	result.m[1][0] = 0;
	result.m[1][1] = 1;
	result.m[1][2] = 0;
	result.m[1][3] = 0;
	result.m[2][0] = 0;
	result.m[2][1] = 0;
	result.m[2][2] = 1;
	result.m[2][3] = 0;
	result.m[3][0] = 0;
	result.m[3][1] = 0;
	result.m[3][2] = 0;
	result.m[3][3] = 1;

	return result;

}

//平行移動行列
Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {

	Matrix4x4 result = {};
	result.m[0][0] = 1.0f;
	result.m[0][1] = 0.0f;
	result.m[0][2] = 0.0f;
	result.m[0][3] = 0.0f;

	result.m[1][0] = 0.0f;
	result.m[1][1] = 1.0f;
	result.m[1][2] = 0.0f;
	result.m[1][3] = 0.0f;

	result.m[2][0] = 0.0f;
	result.m[2][1] = 0.0f;
	result.m[2][2] = 1.0f;
	result.m[2][3] = 0.0f;

	result.m[3][2] = translate.z;
	result.m[3][1] = translate.y;
	result.m[3][0] = translate.x;
	result.m[3][3] = 1.0f;

	return result;
}

//拡大縮小行列
Matrix4x4 MakeScaleMatrix(const Vector3& scale) {
	Matrix4x4 result = {};
	result.m[0][0] = scale.x;
	result.m[0][1] = 0.0f;
	result.m[0][2] = 0.0f;
	result.m[0][3] = 0.0f;

	result.m[1][0] = 0.0f;
	result.m[1][1] = scale.y;
	result.m[1][2] = 0.0f;
	result.m[1][3] = 0.0f;

	result.m[2][0] = 0.0f;
	result.m[2][1] = 0.0f;
	result.m[2][2] = scale.z;
	result.m[2][3] = 0.0f;

	result.m[3][0] = 0.0f;
	result.m[3][1] = 0.0f;
	result.m[3][2] = 0.0f;
	result.m[3][3] = 1.0f;

	return result;
}


//回転行列
Matrix4x4 RotationX(float radian) {
	float cosA = cosf(radian);
	float sinA = sinf(radian);

	Matrix4x4 result = {};
	result.m[0][0] = 1;
	result.m[1][1] = cosA;
	result.m[1][2] = sinA;
	result.m[2][1] = -sinA;
	result.m[2][2] = cosA;
	result.m[3][3] = 1;


	return result;
}


Matrix4x4 RotationY(float radian) {
	float cosA = cosf(radian);
	float sinA = sinf(radian);

	Matrix4x4 result = {};
	result.m[0][0] = cosA;
	result.m[0][2] = -sinA;
	result.m[1][1] = 1;
	result.m[2][0] = sinA;
	result.m[2][2] = cosA;
	result.m[3][3] = 1;


	return result;
}

Matrix4x4 RotationZ(float radian) {
	float cosA = cosf(radian);
	float sinA = sinf(radian);

	Matrix4x4 result = {};
	result.m[0][0] = cosA;
	result.m[0][1] = sinA;
	result.m[1][0] = -sinA;
	result.m[1][1] = cosA;
	result.m[2][2] = 1;
	result.m[3][3] = 1;

	return result;
}

//座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	Vector3 result = {};

	result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0];
	result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1];
	result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];

	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
}

//アフィン変換
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {

	Matrix4x4 makeRotateXMatrix = RotationX(rotate.x);
	Matrix4x4 makeRotateYMatrix = RotationY(rotate.y);
	Matrix4x4 makeRotateZMatrix = RotationZ(rotate.z);

	Matrix4x4 makeRotate = Multiply(Multiply(makeRotateXMatrix, makeRotateYMatrix), makeRotateZMatrix);

	Matrix4x4 matScale = MakeScaleMatrix(scale);
	Matrix4x4 matTranslate = MakeTranslateMatrix(translate);

	return Multiply(Multiply(matScale, makeRotate), matTranslate);
}

//透視投影
Matrix4x4 MakePersectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result = {};


	result.m[0][0] = (1.0f / aspectRatio) * 1.0f / tanf(fovY / 2.0f);
	result.m[1][1] = 1.0f / tanf(fovY / 2.0f);
	result.m[2][2] = (farClip + nearClip) / (farClip - nearClip);
	result.m[2][3] = 1.0f;
	result.m[3][2] = (farClip * nearClip) / -(farClip - nearClip);

	return result;
}

//正射影行列
Matrix4x4 MakeOrthograpicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip) {

	Matrix4x4 result = {};

	result.m[0][0] = 2.0f / (right - left);
	result.m[1][1] = 2.0f / (top - bottom);
	result.m[2][2] = 1.0f / (farClip - nearClip);
	result.m[3][0] = (left + right) / (left - right);
	result.m[3][1] = (top + bottom) / (bottom - top);
	result.m[3][2] = nearClip / (nearClip - farClip);
	result.m[3][3] = 1.0f;
	return result;
}

//ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {

	Matrix4x4 result = {};


	result.m[0][0] = width / 2.0f;
	result.m[1][1] = -(height / 2.0f);
	result.m[2][2] = (maxDepth - minDepth);
	result.m[3][0] = left + (width / 2.0f);
	result.m[3][1] = top + (height / 2.0f);
	result.m[3][2] = minDepth;
	result.m[3][3] = 1.0f;

	return result;
}

//クロス積
Vector3 Cross(const Vector3& v1, const Vector3& v2) {

	Vector3 result;
	result.x = v1.y * v2.z - v1.z * v2.y;
	result.y = v1.z * v2.x - v1.x * v2.z;
	result.z = v1.x * v2.y - v1.y * v2.x;
	return result;

}


void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;
	const uint32_t kSubdivision = 10;
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		float x = -kGridHalfWidth + (xIndex * kGridEvery);
		Vector3 start{ x, 0.0f, -kGridHalfWidth };
		Vector3 end{ x, 0.0f, kGridHalfWidth };
		Vector3 startScreen = Transform(Transform(start, viewProjectionMatrix), viewportMatrix);
		Vector3 endScreen = Transform(Transform(end, viewProjectionMatrix), viewportMatrix);
		Novice::DrawLine(
			int(startScreen.x), int(startScreen.y), int(endScreen.x), int(endScreen.y),
			x == 0.0f ? BLACK : 0xAAAAAAFF);
	}
	for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
		float z = -kGridHalfWidth + (zIndex * kGridEvery);
		Vector3 start{ -kGridHalfWidth, 0.0f, z };
		Vector3 end{ kGridHalfWidth, 0.0f, z };
		Vector3 startScreen = Transform(Transform(start, viewProjectionMatrix), viewportMatrix);
		Vector3 endScreen = Transform(Transform(end, viewProjectionMatrix), viewportMatrix);
		Novice::DrawLine(
			int(startScreen.x), int(startScreen.y), int(endScreen.x), int(endScreen.y),
			z == 0.0f ? BLACK : 0xAAAAAAFF);
	}
}






void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	float pi = std::numbers::pi_v<float>;
	const uint32_t kSubdivision = 12;
	// 経度分割1つ分の角度
	const float kLonEvery = pi * 2.0f / float(kSubdivision);
	// 緯度分割1つ分の角度
	const float kLatEvery = pi / float(kSubdivision);
	// 緯度の方向に分割
	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
		float lat = -pi / 2.0f + kLatEvery * latIndex;
		// 経度の方向に分割しながら線を描く
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery;
			Vector3 a = {
			  sphere.center.x + sphere.radius * cosf(lat) * cosf(lon),
			  sphere.center.y + sphere.radius * sinf(lat),
			  sphere.center.z + sphere.radius * cosf(lat) * sinf(lon) };
			Vector3 b = {
			  sphere.center.x + sphere.radius * cosf(lat + kLatEvery) * cosf(lon),
			  sphere.center.y + sphere.radius * sinf(lat + kLatEvery),
			  sphere.center.z + sphere.radius * cosf(lat + kLatEvery) * sinf(lon) };
			Vector3 c = {
			  sphere.center.x + sphere.radius * cosf(lat) * cosf(lon + kLonEvery),
			  sphere.center.y + sphere.radius * sinf(lat),
			  sphere.center.z + sphere.radius * cosf(lat) * sinf(lon + kLonEvery) };
			// 線を描く
			Vector3 screenA = Transform(Transform(a, viewProjectionMatrix), viewportMatrix);
			Vector3 screenB = Transform(Transform(b, viewProjectionMatrix), viewportMatrix);
			Vector3 screenC = Transform(Transform(c, viewProjectionMatrix), viewportMatrix);
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenB.x), int(screenB.y), color);
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenC.x), int(screenC.y), color);
		}
	}
}

//正射影ベクトル
Vector3 Project(const Vector3& v1, const Vector3& v2) {
	Vector3 result;

	result.x = Dot(v1, Normalize(v2)) * Normalize(v2).x;
	result.y = Dot(v1, Normalize(v2)) * Normalize(v2).y;
	result.z = Dot(v1, Normalize(v2)) * Normalize(v2).z;

	return result;
}
//最近接点
Vector3 ClosestPoint(const Vector3& point, const Segment segment) {


	return Add(Project(Subtract(point, segment.origin), segment.diff), segment.origin);


}

Matrix4x4 MakeViewProjectionMatrix(Vector3 scale, Vector3 rotate, Vector3 translate, Vector3 cameraScale, Vector3 cameraRotate, Vector3 cameraTranslate) {
	Matrix4x4 worldMatrix = MakeAffineMatrix(scale, rotate, translate);
	Matrix4x4 cameraMatrix = MakeAffineMatrix(cameraScale, cameraRotate, cameraTranslate);
	Matrix4x4 viewMatrix = Invers(cameraMatrix);
	Matrix4x4 projectionMatrix = MakePersectiveFovMatrix(0.45f, 720.0f / 1280.0f, 0.1f, 100.0f);
	MatrixScreenPrintf(0, 0, viewMatrix, "viewprojection");
	return (Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix)));

}
//球と球の衝突

bool IsCollision(const Sphere& s1, const Sphere& s2) {
	float distance = Length(Subtract(s2.center, s1.center));
	if (distance <= s1.radius + s2.radius) {

		return true;

	}
	return false;
}

//平面
bool IsCollision(const Sphere& sphere, const Plane& plane) {

	if (sphere.radius >= fabsf(Dot(plane.normal, sphere.center) - plane.distance)) {

		return true;

	}
	return false;
}

//線と平面
bool IsCollision(const Segment& segment, const Plane& plane) {
	float dot = Dot(plane.normal, segment.diff);
	if (dot == 0.0f) {
		return false;
	}
	float t = (plane.distance - Dot(segment.origin, plane.normal)) / dot;

	return t >= 0.0f && t <= 1.0f;
}
//三角形と線
bool IsCollision(const Triangle& triangle, const Segment& segmrnt) {

	Plane plane;
	Vector3 v01 = triangle.vertices[1] - triangle.vertices[0];
	Vector3 v12 = triangle.vertices[2] - triangle.vertices[1];
	Vector3 v20 = triangle.vertices[0] - triangle.vertices[2];

	plane.normal = Normalize(Cross(v01, v12));
	plane.distance = Dot(triangle.vertices[0], plane.normal);

	float t = (plane.distance - Dot(segmrnt.origin, plane.normal)) / Dot(plane.normal, segmrnt.diff);

	Vector3 tb = Multiply(t, segmrnt.diff);

	Vector3 p = segmrnt.origin + tb;

	Vector3 v0p = p - triangle.vertices[0];
	Vector3 v1p = p - triangle.vertices[1];
	Vector3 v2p = p - triangle.vertices[2];

	Vector3 cross01 = Cross(v01, v1p);
	Vector3 cross12 = Cross(v12, v2p);
	Vector3 cross20 = Cross(v20, v0p);

	if (Dot(cross01, plane.normal) >= 0.0f &&
		Dot(cross12, plane.normal) >= 0.0f &&
		Dot(cross20, plane.normal) >= 0.0f) {
		return true;
	}
	return false;
}
//aabbとaabb
bool IsCollision(const AABB& aabb1, const AABB& aabb2) {

	if ((aabb1.min.x <= aabb2.max.x && aabb1.max.x >= aabb2.min.x) &&
		(aabb1.min.y <= aabb2.max.y && aabb1.max.y >= aabb2.min.y) &&
		(aabb1.min.z <= aabb2.max.z && aabb1.max.z >= aabb2.min.z)) {
		return true;
	}

	return false;

}
//aabbと球
bool IsCollision(const AABB& aabb, const Sphere& sphere) {

	Vector3 closestPoint{ std::clamp(sphere.center.x, aabb.min.x, aabb.max.x),
						std::clamp(sphere.center.y, aabb.min.y, aabb.max.y),
						std::clamp(sphere.center.z, aabb.min.z, aabb.max.z) };
	float distance = Length(closestPoint - sphere.center);
	if (distance <= sphere.radius) {
		return true;
	}

	return false;
}
//aabbと線
bool IsCollision(const AABB& aabb, const Segment& segment) {
	float tXMin = (aabb.min.x - segment.origin.x) / segment.diff.x;
	float tYMin = (aabb.min.y - segment.origin.y) / segment.diff.y;
	float tZMin = (aabb.min.z - segment.origin.z) / segment.diff.z;
	float tXMax = (aabb.max.x - segment.origin.x) / segment.diff.x;
	float tYMax = (aabb.max.y - segment.origin.y) / segment.diff.y;
	float tZMax = (aabb.max.z - segment.origin.z) / segment.diff.z;

	float tNearX = min(tXMin, tXMax);
	float tNearY = min(tYMin, tYMax);
	float tNearZ = min(tZMin, tZMax);
	float tFarX = max(tXMin, tXMax);
	float tFarY = max(tYMin, tYMax);
	float tFarZ = max(tZMin, tZMax);


	float tMin = max(max(tNearX, tNearY), tNearZ);

	float tMax = min(min(tFarX, tFarY), tFarZ);

	if (tMin <= tMax) {
		return true;
	}

	return false;
}

Vector3 Perpendicular(const Vector3& vector) {
	if (vector.x != 0.0f || vector.y != 0.0f) {
		return { -vector.y, vector.x, 0.0f };
	}
	return { 0.0f,-vector.z,vector.y };
}

void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {

	Vector3 center = Multiply(plane.distance, plane.normal);
	Vector3 perpendistance[4];

	perpendistance[0] = Normalize(Perpendicular(plane.normal));
	perpendistance[1] = { -perpendistance[0].x,-perpendistance[0].y,-perpendistance[0].z };
	perpendistance[2] = Cross(plane.normal, perpendistance[0]);
	perpendistance[3] = { -perpendistance[2].x,-perpendistance[2].y,perpendistance[2].z };

	Vector3 points[4];
	for (int32_t index = 0; index < 4; index++)
	{
		Vector3 extend = Multiply(2.0f, perpendistance[index]);
		Vector3 point = Add(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}
	Novice::DrawLine(int(points[0].x), int(points[0].y), int(points[2].x), int(points[2].y), color);
	Novice::DrawLine(int(points[1].x), int(points[1].y), int(points[2].x), int(points[2].y), color);
	Novice::DrawLine(int(points[3].x), int(points[3].y), int(points[1].x), int(points[1].y), color);
	Novice::DrawLine(int(points[3].x), int(points[3].y), int(points[0].x), int(points[0].y), color);
}
void DrawTriangle(const Triangle& triangle, const Matrix4x4& viewProjectionMattrix, const Matrix4x4& viewportMatrix, uint32_t color) {

	Vector3 start[3];
	Vector3 end[3];

	for (uint32_t i = 0; i < 3; i++) {
		start[i] = Transform(Transform(triangle.vertices[i], viewProjectionMattrix), viewportMatrix);
		end[i] = Transform(Transform(triangle.vertices[(i + 1) % 3], viewProjectionMattrix), viewportMatrix);
		Novice::DrawLine(int(start[i].x), int(start[i].y), int(end[i].x), int(end[i].y), color);
	}
}
void DrawAABB(const AABB& aabb, const Matrix4x4& viewProjectionMatrix, const Matrix4x4 viewportMatrix, uint32_t color) {
	Vector3 end[12]{};
	Vector3 start[12]{};

	start[0] = Transform(Transform(aabb.min, viewProjectionMatrix), viewportMatrix);
	end[0] = Transform(Transform({ aabb.min.x, aabb.min.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	start[1] = Transform(Transform(aabb.min, viewProjectionMatrix), viewportMatrix);
	end[1] = Transform(Transform({ aabb.min.x, aabb.max.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	start[2] = Transform(Transform({ aabb.min.x, aabb.max.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	end[2] = Transform(Transform({ aabb.min.x, aabb.max.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	start[3] = Transform(Transform({ aabb.min.x, aabb.max.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	end[3] = Transform(Transform({ aabb.min.x, aabb.min.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	start[4] = Transform(Transform({ aabb.min.x, aabb.min.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	end[4] = Transform(Transform({ aabb.max.x, aabb.min.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	start[5] = Transform(Transform(aabb.min, viewProjectionMatrix), viewportMatrix);
	end[5] = Transform(Transform({ aabb.max.x, aabb.min.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	start[6] = Transform(Transform({ aabb.min.x, aabb.max.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	end[6] = Transform(Transform({ aabb.max.x, aabb.max.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	start[7] = Transform(Transform({ aabb.min.x, aabb.max.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	end[7] = Transform(Transform({ aabb.max.x, aabb.max.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	start[8] = Transform(Transform(aabb.max, viewProjectionMatrix), viewportMatrix);
	end[8] = Transform(Transform({ aabb.max.x, aabb.max.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	start[9] = Transform(Transform({ aabb.max.x, aabb.max.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	end[9] = Transform(Transform({ aabb.max.x, aabb.min.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	start[10] = Transform(Transform({ aabb.max.x, aabb.min.y, aabb.min.z }, viewProjectionMatrix), viewportMatrix);
	end[10] = Transform(Transform({ aabb.max.x, aabb.min.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	start[11] = Transform(Transform({ aabb.max.x, aabb.min.y, aabb.max.z }, viewProjectionMatrix), viewportMatrix);
	end[11] = Transform(Transform(aabb.max, viewProjectionMatrix), viewportMatrix);

	for (uint32_t i = 0; i < 12; i++) {
		Novice::DrawLine(int(start[i].x), int(start[i].y), int(end[i].x), int(end[i].y), color);
	}
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };


	Vector3 cameraTranslate{ 0.0f,0.0f,-10.0f };
	Vector3 cameraRotato{ 0.0f,0.0f,0.0f };

	Segment segmrnt;
	segmrnt.diff = { 0.45f,0.78f,0.0f };
	segmrnt.origin = { 1.0f,0.58f,0.0f };


	AABB aabb{
				.min{-0.5f, -0.5f, -0.5f},
				.max{0.0f, 0.0f, 0.0f},
	};
	Vector3 rotate{};
	Vector3 translate{};
	Vector3 Scale = { 1.0f,1.0f,1.0f };

	Vector3 start{};
	Vector3 end{};

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		Matrix4x4 worldMatrix = MakeAffineMatrix(Scale, rotate, translate);
		Matrix4x4 cameraMatrix = MakeAffineMatrix(Scale, cameraRotato, cameraTranslate);
		Matrix4x4 viewMatrix = Invers(cameraMatrix);
		Matrix4x4 projectionMatrix = MakePersectiveFovMatrix(0.45f, float(kWindowWith) / float(kWindowHigat), 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = Multiply(worldMatrix, Multiply(viewMatrix, projectionMatrix));
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWith), float(kWindowHigat), 0.0f, 1.0f);



		ImGui::DragFloat3("sphere.center", &segmrnt.diff.x, 0.01f);
		ImGui::DragFloat("sphere.origin", &segmrnt.origin.x, 0.01f);
		ImGui::DragFloat3("aabb.min", &aabb.min.x, 0.01f);
		ImGui::DragFloat3("aabb.max", &aabb.max.x, 0.01f);
		ImGui::DragFloat3("cameraRotate", &cameraRotato.x, 0.01f);
		ImGui::DragFloat3("cameraTranslate", &cameraTranslate.x, 0.01f);

		start = Transform(Transform(segmrnt.diff, viewProjectionMatrix), viewportMatrix);
		end = Transform(Transform(segmrnt.origin, viewProjectionMatrix), viewportMatrix);

		///
		///
		/// ↑更新処理ここまで
		///




		///
		/// ↓描画処理ここから
		///
		DrawGrid(viewProjectionMatrix, viewportMatrix);

		Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, WHITE);
		DrawAABB(aabb, viewProjectionMatrix, viewportMatrix, WHITE);


		if (IsCollision(aabb, segmrnt)) {
			Novice::DrawLine((int)start.x, (int)start.y, (int)end.x, (int)end.y, RED);

		}


		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
