# =============================================================================
# eda_report.R — Sinh hình ảnh và số liệu cho Section 2 (Dữ liệu & EDA)
# Output ảnh: report/img/fig_s2_0*.png
# Chạy từ thư mục gốc project:
#   Rscript src/reproduction/eda_report.R
# =============================================================================

source("src/00_setup.R")

# Ghi đè thư mục đầu ra ảnh về report/img/
IMG_DIR <- file.path(PROJECT_ROOT, "report", "img")
if (!dir.exists(IMG_DIR)) dir.create(IMG_DIR, recursive = TRUE)

img <- function(filename, expr, width = 1000, height = 600, res = 120) {
  path <- file.path(IMG_DIR, filename)
  png(path, width = width, height = height, res = res)
  eval(expr)
  dev.off()
  cat("  [saved]", filename, "\n")
}

# =============================================================================
# 1. Đọc & chuẩn bị dữ liệu
# =============================================================================
csv_path <- file.path(DATA_DIR, "AEP_hourly.csv")
if (!file.exists(csv_path)) stop("Không tìm thấy: ", csv_path)

raw <- read_csv(csv_path, show_col_types = FALSE) %>%
  rename(datetime = Datetime, value = AEP_MW) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "America/New_York")) %>%
  arrange(datetime) %>%
  distinct(datetime, .keep_all = TRUE) %>%
  filter(!is.na(datetime), !is.na(value))

START_DATE <- ymd_hms("2017-06-01 00:00:00")
END_DATE   <- ymd_hms("2017-08-31 23:00:00")

dat <- raw %>% filter(datetime >= START_DATE, datetime <= END_DATE)

# Nội suy giờ thiếu
n_expected <- as.integer(difftime(END_DATE, START_DATE, units = "hours")) + 1
n_missing  <- n_expected - nrow(dat)
full_seq   <- seq(START_DATE, END_DATE, by = "hour")
dat <- data.frame(datetime = full_seq) %>%
  left_join(dat, by = "datetime")
dat$value <- approx(seq_along(dat$value), dat$value, seq_along(dat$value))$y

dat$hour       <- hour(dat$datetime)
dat$wday       <- wday(dat$datetime, label = TRUE, abbr = TRUE, locale = "C")
dat$is_weekend <- as.integer(wday(dat$datetime) %in% c(1, 7))

y <- ts(dat$value, frequency = 24)

# =============================================================================
# 2. Thống kê mô tả
# =============================================================================
cat("\n=== THỐNG KÊ MÔ TẢ ===\n")
cat(sprintf("N        = %d\n",       length(y)))
cat(sprintf("Mean     = %.2f MW\n",  mean(y)))
cat(sprintf("SD       = %.2f MW\n",  sd(y)))
cat(sprintf("Min      = %.2f MW\n",  min(y)))
cat(sprintf("Q1       = %.2f MW\n",  quantile(y, 0.25)))
cat(sprintf("Median   = %.2f MW\n",  median(y)))
cat(sprintf("Q3       = %.2f MW\n",  quantile(y, 0.75)))
cat(sprintf("Max      = %.2f MW\n",  max(y)))
cat(sprintf("Missing  = %d giờ (đã nội suy)\n", n_missing))

# =============================================================================
# 3. Plots
# =============================================================================

# --- Fig 1: Toàn chuỗi -------------------------------------------------------
img("fig_s2_01_full_series.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  plot(dat$datetime, y, type = "l", col = "steelblue", lwd = 0.8,
       main = "Tiêu thụ điện AEP theo giờ (01/06 – 31/08/2017)",
       xlab = "Thời gian", ylab = "Tiêu thụ (MW)",
       cex.main = 1.1, cex.lab = 0.95)
  grid(col = "gray90", lty = 1)
}), width = 1400, height = 500)

# --- Fig 2: Zoom 2 tuần -------------------------------------------------------
img("fig_s2_02_weekly_zoom.png", expression({
  idx <- 1:336    # 14 ngày = 14*24 = 336 giờ
  par(mar = c(4, 4.5, 3, 1))
  plot(dat$datetime[idx], y[idx], type = "l", col = "darkred", lwd = 1.2,
       main = "Zoom 2 tuần đầu — chu kỳ ngày & tuần",
       xlab = "Thời gian", ylab = "Tiêu thụ (MW)",
       cex.main = 1.1, cex.lab = 0.95)
  # Vạch dọc mỗi 24h
  abline(v = dat$datetime[seq(1, 336, 24)], col = "gray60", lty = 2, lwd = 0.7)
  # Tô màu cuối tuần
  for (i in seq(1, 336, 1)) {
    if (dat$is_weekend[i] == 1) {
      rect(dat$datetime[max(1,i)], par("usr")[3],
           dat$datetime[min(336,i+1)], par("usr")[4],
           col = rgb(1, 0.8, 0.8, 0.15), border = NA)
    }
  }
  legend("topright", c("MW theo giờ", "Ranh giới ngày", "Cuối tuần"),
         col = c("darkred", "gray60", rgb(1,0.5,0.5,0.4)),
         lty = c(1,2,NA), pch = c(NA,NA,15), pt.cex = 1.5,
         lwd = c(1.2,0.7,NA), cex = 0.85, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 1200, height = 450)

# --- Fig 3: Boxplot theo giờ -------------------------------------------------
img("fig_s2_03_diurnal_box.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  boxplot(value ~ hour, data = dat,
          main = "Phân bố tiêu thụ điện theo giờ trong ngày",
          xlab = "Giờ trong ngày (0 = 0:00, 23 = 23:00)",
          ylab = "Tiêu thụ (MW)",
          col = "lightblue", border = "steelblue",
          cex.main = 1.1, cex.lab = 0.95,
          names = 0:23, las = 1)
  grid(nx = NA, ny = NULL, col = "gray90", lty = 1)
}), width = 1100, height = 600)

# --- Fig 4: Ngày thường vs cuối tuần (median daily profile) ------------------
img("fig_s2_04_weekday_vs_weekend.png", expression({
  # Tính median MW theo giờ cho mỗi nhóm
  wd_med <- tapply(dat$value[dat$is_weekend == 0],
                   dat$hour[dat$is_weekend == 0], median)
  we_med <- tapply(dat$value[dat$is_weekend == 1],
                   dat$hour[dat$is_weekend == 1], median)
  ylim_r  <- range(c(wd_med, we_med)) + c(-200, 200)
  par(mar = c(4, 4.5, 3, 1))
  plot(0:23, wd_med, type = "b", col = "steelblue", lwd = 2, pch = 16,
       ylim = ylim_r,
       main = "Profile tiêu thụ trung bình: Ngày thường vs Cuối tuần",
       xlab = "Giờ trong ngày", ylab = "Tiêu thụ trung vị (MW)",
       cex.main = 1.1, cex.lab = 0.95, xaxt = "n")
  axis(1, at = 0:23, labels = 0:23)
  lines(0:23, we_med, type = "b", col = "salmon", lwd = 2, pch = 17)
  legend("bottomright", c("Ngày thường (Mon–Fri)", "Cuối tuần (Sat–Sun)"),
         col = c("steelblue","salmon"), lty = 1, lwd = 2, pch = c(16,17),
         cex = 0.9, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 900, height = 600)

# --- Fig 5: Decompose cộng tính ----------------------------------------------
img("fig_s2_05_decompose.png", expression({
  plot(decompose(y, type = "additive"))
}), width = 1000, height = 800)

# --- Fig 6: ACF / PACF chuỗi gốc --------------------------------------------
img("fig_s2_06_acf_pacf.png", expression({
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  acf(y,  lag.max = 72,
      main = "ACF — chuỗi gốc (lag.max = 72)",
      cex.main = 1.0, cex.lab = 0.9)
  pacf(y, lag.max = 72,
       main = "PACF — chuỗi gốc (lag.max = 72)",
       cex.main = 1.0, cex.lab = 0.9)
  par(mfrow = c(1, 1))
}), width = 1200, height = 600)

# --- Fig 7: Sau sai phân mùa lag=24 -----------------------------------------
y_sd <- diff(y, lag = 24)
img("fig_s2_07_seasonal_diff.png", expression({
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  plot(y_sd, col = "darkgreen", lwd = 0.8,
       main = "Chuỗi sau sai phân mùa (lag = 24)",
       ylab = "ΔMW", xlab = "Thời gian (đơn vị 24h)",
       cex.main = 1.0)
  abline(h = 0, lty = 2, col = "red")
  grid(col = "gray90", lty = 1)
  acf(y_sd, lag.max = 72,
      main = "ACF sau sai phân mùa",
      cex.main = 1.0)
  par(mfrow = c(1, 1))
}), width = 1200, height = 700)

# =============================================================================
# 4. Kiểm định tính dừng
# =============================================================================
cat("\n=== ADF TEST (chuỗi gốc) ===\n")
adf_orig <- adf.test(y)
print(adf_orig)
if (adf_orig$p.value < 0.05) {
  cat("Kết luận: Bác bỏ H0 (không dừng) => chuỗi DỪNG (p <", round(adf_orig$p.value,4), ")\n")
} else {
  cat("Kết luận: Không bác bỏ H0 => chuỗi KHÔNG DỪNG (p =", round(adf_orig$p.value,4), ")\n")
}

cat("\n=== KPSS TEST (chuỗi gốc) ===\n")
kpss_orig <- kpss.test(y)
print(kpss_orig)
if (kpss_orig$p.value < 0.05) {
  cat("Kết luận: Bác bỏ H0 (dừng) => chuỗi KHÔNG DỪNG (p <", round(kpss_orig$p.value,4), ")\n")
} else {
  cat("Kết luận: Không bác bỏ H0 => chuỗi DỪNG (p =", round(kpss_orig$p.value,4), ")\n")
}

cat("\n=== ADF TEST (sau sai phân mùa lag=24) ===\n")
adf_diff <- adf.test(y_sd)
print(adf_diff)
if (adf_diff$p.value < 0.05) {
  cat("Kết luận: Bác bỏ H0 => chuỗi sai phân mùa DỪNG (p <", round(adf_diff$p.value,4), ")\n")
} else {
  cat("Kết luận: Không bác bỏ H0 => chuỗi sai phân mùa KHÔNG DỪNG (p =", round(adf_diff$p.value,4), ")\n")
}

cat("\n✓ Hoàn tất eda_report.R — ảnh đã lưu vào report/img/\n")
