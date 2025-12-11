#include <gtest/gtest.h>
#include "TWDeckWfmFilter.h"
#include "TWDeckWfm.h"
#include <cmath>

// Test Plan:
// 1. Normal cases for each filter shape (kSqr, kTrngl, kGauss, kDiff)
// 2. Edge cases: N <= bw, N <= 2*bw, N <= 8*bw, bw = 0, negative N
// 3. Error handling: unknown filter shape
// 4. Check correct size, bandwidth, shift, and waveform values

TEST(TWDeckWfmFilterTest, SquareNormal) {
  int N = 10;
  double bw = 5;
  TWDeckWfmFilter filter(N, wdeck::kSqr, bw);
  EXPECT_EQ(filter.GetBandwidth(), bw);
  EXPECT_EQ(filter.GetShift(), 0.5 * bw);
  EXPECT_EQ(filter.GetSize(), N);
  for (int i = 0; i < bw; ++i) {
    EXPECT_DOUBLE_EQ(filter.GetWfm().at(i), 1.0 / bw);
  }
}

TEST(TWDeckWfmFilterTest, SquareEdgeNLessThanBw) {
  int N = 4;
  double bw = 5;
  TWDeckWfmFilter filter(N, wdeck::kSqr, bw);
  EXPECT_EQ(filter.GetSize(), 2 * bw);
  for (int i = 0; i < bw; ++i) {
    EXPECT_DOUBLE_EQ(filter.GetWfm().at(i), 1.0 / bw);
  }
}

TEST(TWDeckWfmFilterTest, TriangleNormal) {
  int N = 20;
  double bw = 5;
  TWDeckWfmFilter filter(N, wdeck::kTrngl, bw);
  EXPECT_EQ(filter.GetBandwidth(), bw);
  EXPECT_EQ(filter.GetShift(), bw);
  EXPECT_EQ(filter.GetSize(), N);
  double h = 1.0 / bw;
  double slope = h / bw;
  for (int i = 0; i < 2 * bw; ++i) {
    double expected = (i <= bw) ? i * slope : h - slope * (i - bw);
    EXPECT_NEAR(filter.GetWfm().at(i), expected, 1e-12);
  }
}

TEST(TWDeckWfmFilterTest, TriangleEdgeNLessThanBw) {
  int N = 4;
  double bw = 5;
  TWDeckWfmFilter filter(N, wdeck::kTrngl, bw);
  EXPECT_EQ(filter.GetSize(), 3 * bw);
}

TEST(TWDeckWfmFilterTest, GaussNormal) {
  int N = 50;
  double bw = 5;
  TWDeckWfmFilter filter(N, wdeck::kGauss, bw);
  EXPECT_EQ(filter.GetBandwidth(), bw);
  EXPECT_EQ(filter.GetShift(), 4 * bw);
  EXPECT_EQ(filter.GetSize(), N);
  for (int i = 0; i < 8 * bw; ++i) {
    double expected = TMath::Gaus(i, 4 * bw, bw, true);
    EXPECT_NEAR(filter.GetWfm().at(i), expected, 1e-12);
  }
}

TEST(TWDeckWfmFilterTest, GaussEdgeNLessThan8Bw) {
  int N = 10;
  double bw = 2;
  TWDeckWfmFilter filter(N, wdeck::kGauss, bw);
  EXPECT_EQ(filter.GetSize(), 8 * bw);
}

TEST(TWDeckWfmFilterTest, DiffNormal) {
  int N = 20;
  double bw = 5;
  TWDeckWfmFilter filter(N, wdeck::kDiff, bw);
  EXPECT_EQ(filter.GetBandwidth(), bw);
  EXPECT_EQ(filter.GetShift(), bw);
  EXPECT_EQ(filter.GetSize(), N);
  for (int i = 0; i < 2 * bw; ++i) {
    double expected = (i < bw) ? 1.0 / bw : -1.0 / bw;
    EXPECT_DOUBLE_EQ(filter.GetWfm().at(i), expected);
  }
}

TEST(TWDeckWfmFilterTest, DiffEdgeNLessThan2Bw) {
  int N = 8;
  double bw = 5;
  TWDeckWfmFilter filter(N, wdeck::kDiff, bw);
  EXPECT_EQ(filter.GetSize(), 3 * bw);
}

TEST(TWDeckWfmFilterTest, UnknownShape) {
  int N = 10;
  double bw = 5;
  TWDeckWfmFilter filter(N, static_cast<wdeck::EFltrShape>(999), bw);
  // Should print error, but no exception thrown. Optionally, check for unchanged state.
}

TEST(TWDeckWfmFilterTest, ZeroBandwidth) {
  int N = 10;
  double bw = 0;
  TWDeckWfmFilter filter(N, wdeck::kSqr, bw);
  EXPECT_EQ(filter.GetBandwidth(), 0);
  EXPECT_EQ(filter.GetShift(), 0);
  // No waveform values should be set
}

