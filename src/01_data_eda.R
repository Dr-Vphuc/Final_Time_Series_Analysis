# =============================================================================
# 01_data_eda.R — Tải dữ liệu, EDA, kiểm định tính dừng
# Tương ứng NGÀY 1 trong timeline.
#
# DỮ LIỆU: PJM Hourly Energy Consumption (Kaggle)
#   https://www.kaggle.com/datasets/robikscube/hourly-energy-consumption
#
# Tải file AEP_hourly.csv và đặt vào: data/AEP_hourly.csv
# =============================================================================

source("src/00_setup.R")

# --- 1. Đọc dữ liệu ----------------------------------------------------------
csv_path <- file.path(DATA_DIR, "AEP_hourly.csv")

if (!file.exists(csv_path)) {
  stop("\n[!] Chưa tìm thấy ", csv_path,
       "\n    Hãy tải AEP_hourly.csv từ Kaggle và đặt vào thư mục data/.")
}

raw <- read_csv(csv_path, show_col_types = FALSE)
raw <- raw %>%
  rename(datetime = Datetime, value = AEP_MW) %>%
  mutate(datetime = ymd_hms(datetime)) %>%
  arrange(datetime) %>%
  distinct(datetime, .keep_all = TRUE)

cat("Dữ liệu thô:", nrow(raw), "quan sát từ",
    as.character(min(raw$datetime)), "đến",
    as.character(max(raw$datetime)), "\n")

# --- 2. Cắt phạm vi 3 tháng mùa hè -------------------------------------------
# Chọn 06/2017 - 08/2017 (~2.184 quan sát)
START_DATE <- ymd_hms("2017-06-01 00:00:00")
END_DATE   <- ymd_hms("2017-08-31 23:00:00")

dat <- raw %>%
  filter(datetime >= START_DATE, datetime <= END_DATE)

cat("Sau cắt phạm vi:", nrow(dat), "quan sát\n")

# Kiểm missing
n_missing_hours <- as.numeric(difftime(END_DATE, START_DATE, units = "hours")) + 1 - nrow(dat)
cat("Số giờ thiếu:", n_missing_hours, "\n")

# Nội suy nếu thiếu (đơn giản — linear)
full_seq <- seq(START_DATE, END_DATE, by = "hour")
dat <- data.frame(datetime = full_seq) %>%
  left_join(dat, by = "datetime")
dat$value <- approx(seq_along(dat$value), dat$value, seq_along(dat$value))$y

# --- 3. Tạo ts object --------------------------------------------------------
y <- ts(dat$value, frequency = 24)   # frequency = 24 vì chu kỳ ngày
cat("Chuỗi y có", length(y), "quan sát, tần số =", frequency(y), "\n\n")

# --- 4. EDA ------------------------------------------------------------------

# (a) Plot toàn chuỗi
save_plot("01a_full_series.png", expression({
  plot(y, main = "Tiêu thụ điện AEP (06-08/2017)",
       ylab = "MW", xlab = "Ngày (đơn vị 24 giờ)", col = "steelblue")
}))

# (b) Zoom 1 tuần
save_plot("01b_one_week.png", expression({
  plot(y[1:168], type = "l",
       main = "Zoom 1 tuần đầu — thấy chu kỳ ngày + tuần",
       ylab = "MW", xlab = "Giờ", col = "darkred")
  abline(v = seq(0, 168, 24), col = "gray", lty = 2)
}))

# (c) Boxplot theo giờ trong ngày
dat$hour <- hour(dat$datetime)
save_plot("01c_boxplot_hour.png", expression({
  boxplot(value ~ hour, data = dat,
          main = "Phân bố tiêu thụ theo giờ trong ngày",
          xlab = "Giờ", ylab = "MW", col = "lightblue")
}))

# (d) Boxplot weekday vs weekend
dat$is_weekend <- as.integer(wday(dat$datetime) %in% c(1, 7))
save_plot("01d_boxplot_weekend.png", expression({
  boxplot(value ~ is_weekend, data = dat,
          main = "Tiêu thụ: Ngày thường (0) vs Cuối tuần (1)",
          xlab = "is_weekend", ylab = "MW", col = c("lightgreen", "salmon"))
}))

# (e) Decompose
save_plot("01e_decompose_additive.png", expression({
  plot(decompose(y, type = "additive"))
}))

save_plot("01f_decompose_multiplicative.png", expression({
  plot(decompose(y, type = "multiplicative"))
}))

# (f) ACF / PACF với lag.max = 72 (3 ngày)
save_plot("01g_acf_pacf.png", expression({
  par(mfrow = c(1, 2))
  acf(y,  lag.max = 72, main = "ACF (lag = 72)")
  pacf(y, lag.max = 72, main = "PACF (lag = 72)")
  par(mfrow = c(1, 1))
}))

# --- 5. Kiểm định tính dừng --------------------------------------------------
cat("\n=== KIỂM ĐỊNH TÍNH DỪNG ===\n")

cat("\n[ADF test]\n")
print(adf.test(y))

cat("\n[KPSS test]\n")
print(kpss.test(y))

# Sai phân mùa
y_diff <- diff(y, lag = 24)
cat("\n[ADF sau sai phân mùa lag=24]\n")
print(adf.test(y_diff))

save_plot("01h_diff_seasonal.png", expression({
  par(mfrow = c(2, 1))
  plot(y_diff, main = "Sau sai phân mùa lag=24",
       ylab = "MW", col = "darkgreen")
  acf(y_diff, lag.max = 72, main = "ACF của chuỗi sai phân mùa")
  par(mfrow = c(1, 1))
}))

# --- 6. Lưu dữ liệu đã xử lý -------------------------------------------------
saveRDS(list(y = y, dat = dat),
        file = file.path(DATA_DIR, "processed_data.rds"))

cat("\n✓ Hoàn tất Ngày 1. Dữ liệu đã lưu vào processed_data.rds\n")
cat("  Tiếp theo: chạy 02_model1_indicator.R\n")
