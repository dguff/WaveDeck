#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include "TGraph.h"
#include "TWDeckUtils.h"
#include "TF1.h"

// linspace tests
TEST(LinspaceTest, NormalCase) {
    auto result = linspace<double>(0.0, 1.0, 5);
    std::vector<double> expected = {0.0, 0.25, 0.5, 0.75, 1.0};
    ASSERT_EQ(result.size(), expected.size());
    for (size_t i = 0; i < result.size(); ++i) {
        EXPECT_NEAR(result[i], expected[i], 1e-9);
    }
}

TEST(LinspaceTest, SinglePoint) {
    auto result = linspace<double>(2.0, 2.0, 1);
    ASSERT_EQ(result.size(), 1);
    EXPECT_DOUBLE_EQ(result[0], 2.0);
}

TEST(LinspaceTest, TwoPoints) {
    auto result = linspace<double>(-1.0, 1.0, 2);
    ASSERT_EQ(result.size(), 2);
    EXPECT_DOUBLE_EQ(result[0], -1.0);
    EXPECT_DOUBLE_EQ(result[1], 1.0);
}

// g_integral tests
TEST(GIntegralTest, LinearGraph) {
    double x[] = {0, 1, 2};
    double y[] = {0, 1, 2};
    TGraph g(3, x, y);
    double result = g_integral(&g, 0, 2);
    // Integral of y=x from 0 to 2 is 2
    EXPECT_NEAR(result, 2.0, 0.1);
}

TEST(GIntegralTest, ConstantGraph) {
    double x[] = {0, 1, 2};
    double y[] = {5, 5, 5};
    TGraph g(3, x, y);
    double result = g_integral(&g, 0, 2);
    // Integral of y=5 from 0 to 2 is 10
    EXPECT_NEAR(result, 10.0, 0.1);
}

// g_scale_Y tests
TEST(GScaleYTest, ScaleUp) {
    double x[] = {0, 1, 2};
    double y[] = {1, 2, 3};
    TGraph g(3, x, y);
    g_scale_Y(&g, 2.0);
    for (int i = 0; i < g.GetN(); ++i) {
        EXPECT_DOUBLE_EQ(g.GetY()[i], y[i] * 2.0);
    }
}

TEST(GScaleYTest, ScaleDown) {
    double x[] = {0, 1, 2};
    double y[] = {2, 4, 6};
    TGraph g(3, x, y);
    g_scale_Y(&g, 0.5);
    for (int i = 0; i < g.GetN(); ++i) {
        EXPECT_DOUBLE_EQ(g.GetY()[i], y[i] * 0.5);
    }
}

// g_scale_X tests
TEST(GScaleXTest, ScaleUp) {
    double x[] = {1, 2, 3};
    double y[] = {0, 0, 0};
    TGraph g(3, x, y);
    g_scale_X(&g, 3.0);
    for (int i = 0; i < g.GetN(); ++i) {
        EXPECT_DOUBLE_EQ(g.GetX()[i], x[i] * 3.0);
    }
}

TEST(GScaleXTest, ScaleDown) {
    double x[] = {2, 4, 6};
    double y[] = {0, 0, 0};
    TGraph g(3, x, y);
    g_scale_X(&g, 0.5);
    for (int i = 0; i < g.GetN(); ++i) {
        EXPECT_DOUBLE_EQ(g.GetX()[i], x[i] * 0.5);
    }
}

// g_find_x tests
TEST(GFindXTest, FindXLinear) {
    double x[] = {0, 1, 2};
    double y[] = {0, 1, 2};
    TGraph g(3, x, y);
    double result = g_find_x(&g, 1.0, 0, 2);
    EXPECT_NEAR(result, 1.0, 0.05);
}

TEST(GFindXTest, FindXConstant) {
    double x[] = {0, 1, 2};
    double y[] = {5, 5, 5};
    TGraph g(3, x, y);
    double result = g_find_x(&g, 5.0, 0, 2);
    // Any x in [0,2] is valid, but function should return something in range
    EXPECT_GE(result, 0.0);
    EXPECT_LE(result, 2.0);
}

// Edge case: N = 0 for linspace (should handle gracefully)
TEST(LinspaceTest, ZeroPoints) {
    auto result = linspace<double>(0.0, 1.0, 0);
    EXPECT_EQ(result.size(), 0);
}

