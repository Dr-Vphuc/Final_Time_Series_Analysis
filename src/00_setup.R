# =============================================================================
# 00_setup.R — Cài đặt môi trường & hàm tiện ích dùng chung
# Chạy 1 lần đầu tiên hoặc khi mở project mới.
# =============================================================================

# --- 1. Cài package nếu chưa có ----------------------------------------------
required_pkgs <- c(
  "tseries",   # adf.test, kpss.test
  "TSA",       # eacf, harmonic, season, runs
  "astsa",     # sarima, sarima.for
  "forecast",  # auto.arima, fourier, ma, Arima
  "ggplot2",   # vẽ đẹp
  "dplyr",     # xử lý data
  "lubridate", # xử lý date-time
  "readr"      # đọc csv nhanh
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

invisible(lapply(required_pkgs, library, character.only = TRUE))

# --- 2. Cấu hình project -----------------------------------------------------
PROJECT_ROOT <- "f:/Learning/Time_Series_Analysis"
DATA_DIR     <- file.path(PROJECT_ROOT, "data")
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "output")

if (!dir.exists(DATA_DIR))   dir.create(DATA_DIR,   recursive = TRUE)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

set.seed(42)

# --- 3. Hàm tiện ích ---------------------------------------------------------

# Tính các metrics dự báo
compute_metrics <- function(actual, predicted) {
  err  <- actual - predicted
  rmse <- sqrt(mean(err^2, na.rm = TRUE))
  mae  <- mean(abs(err), na.rm = TRUE)
  mape <- mean(abs(err / actual), na.rm = TRUE) * 100
  data.frame(RMSE = rmse, MAE = mae, MAPE = mape)
}

# Lưu plot ra PNG
save_plot <- function(filename, expr, width = 1000, height = 600) {
  png(file.path(OUTPUT_DIR, filename), width = width, height = height, res = 100)
  eval(expr)
  dev.off()
  cat("  → đã lưu:", filename, "\n")
}

# Chia train/test theo tỷ lệ (theo thứ tự thời gian)
train_test_split <- function(y, train_ratio = 0.8) {
  n        <- length(y)
  n_train  <- floor(n * train_ratio)
  list(
    train      = y[1:n_train],
    test       = y[(n_train + 1):n],
    n_train    = n_train,
    n_test     = n - n_train
  )
}

cat("✓ Setup hoàn tất.\n")
cat("  - DATA_DIR:  ", DATA_DIR,   "\n")
cat("  - OUTPUT_DIR:", OUTPUT_DIR, "\n")
