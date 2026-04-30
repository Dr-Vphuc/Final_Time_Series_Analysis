# =============================================================================
# models_report.R — Chạy 3 mô hình, sinh ảnh & số liệu cho Chapter 4
# Output ảnh: report/img/fig_s4_0*.png  và  fig_s5_0*.png
# Chạy từ thư mục gốc project:
#   Rscript src/reproduction/models_report.R
# =============================================================================

source("src/00_setup.R")

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
# 1. Đọc & chuẩn bị dữ liệu (giống eda_report.R)
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

full_seq <- seq(START_DATE, END_DATE, by = "hour")
dat <- data.frame(datetime = full_seq) %>%
  left_join(dat, by = "datetime")
dat$value <- approx(seq_along(dat$value), dat$value, seq_along(dat$value))$y

dat$hour <- hour(dat$datetime)
dat$wday <- wday(dat$datetime)          # 1 = Chủ nhật, 7 = Thứ bảy

# =============================================================================
# 2. Train / Test split  (80 / 20)
# =============================================================================
N       <- nrow(dat)
n_train <- floor(N * 0.8)
n_test  <- N - n_train

dat_train <- dat[1:n_train, ]
dat_test  <- dat[(n_train + 1):N, ]

cat("\n=== TRAIN/TEST SPLIT ===\n")
cat(sprintf("N = %d | Train = %d (%.1f%%) | Test = %d (%.1f%%)\n",
            N, n_train, 100 * n_train / N, n_test, 100 * n_test / N))
cat(sprintf("Train: %s  ->  %s\n",
            format(dat_train$datetime[1],   "%Y-%m-%d %H:%M"),
            format(tail(dat_train$datetime, 1), "%Y-%m-%d %H:%M")))
cat(sprintf("Test:  %s  ->  %s\n",
            format(dat_test$datetime[1],    "%Y-%m-%d %H:%M"),
            format(tail(dat_test$datetime,  1), "%Y-%m-%d %H:%M")))

# Fig 4.1 — vùng train/test
img("fig_s4_01_train_test_split.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  plot(dat$datetime, dat$value, type = "l", col = "gray70", lwd = 0.6,
       main = "Phân chia tập huấn luyện / kiểm tra (80 / 20)",
       xlab = "Thời gian", ylab = "Tiêu thụ (MW)", cex.main = 1.0)
  # tô vùng test
  polygon(c(dat_test$datetime[1], dat_test$datetime,
            tail(dat_test$datetime, 1), dat_test$datetime[1]),
          c(par("usr")[3], dat_test$value,
            par("usr")[3], par("usr")[3]),
          col = rgb(1, 0.6, 0.2, 0.25), border = NA)
  lines(dat_train$datetime, dat_train$value, col = "steelblue", lwd = 0.8)
  lines(dat_test$datetime,  dat_test$value,  col = "darkorange", lwd = 0.8)
  abline(v = dat_test$datetime[1], col = "red", lty = 2, lwd = 1.5)
  legend("topright",
         c(sprintf("Tập huấn luyện (n=%d)", n_train),
           sprintf("Tập kiểm tra   (n=%d)", n_test)),
         col = c("steelblue", "darkorange"), lty = 1, lwd = 1.5,
         cex = 0.85, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 1400, height = 500)

# =============================================================================
# 3. MODEL 1 — Indicator Regression
# =============================================================================
cat("\n=== MODEL 1: INDICATOR REGRESSION ===\n")

dat_train$hour_f <- factor(dat_train$hour)
dat_train$wday_f <- factor(dat_train$wday)
dat_test$hour_f  <- factor(dat_test$hour, levels = levels(dat_train$hour_f))
dat_test$wday_f  <- factor(dat_test$wday, levels = levels(dat_train$wday_f))

t1 <- proc.time()
m1 <- lm(value ~ hour_f + wday_f, data = dat_train)
t1_elapsed <- (proc.time() - t1)["elapsed"]

s1 <- summary(m1)
cat(sprintf("R²       = %.4f\n", s1$r.squared))
cat(sprintf("Adj. R²  = %.4f\n", s1$adj.r.squared))
cat(sprintf("AIC      = %.2f\n", AIC(m1)))
cat(sprintf("Time     = %.3f s\n", t1_elapsed))

pred_m1_train <- fitted(m1)
pred_m1_test  <- predict(m1, newdata = dat_test)
resid_m1      <- residuals(m1)

m1_tr <- compute_metrics(dat_train$value, pred_m1_train)
m1_te <- compute_metrics(dat_test$value,  pred_m1_test)
lb_m1 <- Box.test(resid_m1, lag = 48, type = "Ljung-Box")$p.value

cat(sprintf("Train:  RMSE=%.2f  MAE=%.2f  MAPE=%.2f%%\n",
            m1_tr$RMSE, m1_tr$MAE, m1_tr$MAPE))
cat(sprintf("Test:   RMSE=%.2f  MAE=%.2f  MAPE=%.2f%%\n",
            m1_te$RMSE, m1_te$MAE, m1_te$MAPE))
cat(sprintf("Ljung-Box Q(48) p = %.4f\n", lb_m1))

img("fig_s4_02_model1_fit.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  plot(dat_test$datetime, dat_test$value, type = "l", col = "gray50",
       main = "Mô hình 1 (Indicator): Dự báo vs Thực tế — tập kiểm tra",
       xlab = "Thời gian", ylab = "Tiêu thụ (MW)", cex.main = 1.0)
  lines(dat_test$datetime, pred_m1_test, col = "steelblue", lwd = 1.5)
  legend("topright", c("Thực tế", "Dự báo M1"),
         col = c("gray50", "steelblue"), lty = 1, lwd = c(1, 1.5),
         cex = 0.85, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 1200, height = 450)

img("fig_s4_03_model1_residuals.png", expression({
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  plot(dat_train$datetime, resid_m1, type = "l", col = "steelblue",
       main = "Phần dư M1 theo thời gian",
       xlab = "Thời gian", ylab = "Phần dư (MW)", cex.main = 0.95)
  abline(h = 0, col = "red", lty = 2)
  acf(resid_m1, lag.max = 72, main = "ACF phần dư M1", cex.main = 0.95)
  qqnorm(resid_m1, main = "Q-Q phần dư M1", cex.main = 0.95)
  qqline(resid_m1, col = "red")
  par(mfrow = c(1, 1))
}), width = 1400, height = 500)

# =============================================================================
# 4. MODEL 2 — Fourier Regression  (K = 1 .. 4)
# =============================================================================
cat("\n=== MODEL 2: FOURIER REGRESSION ===\n")

y_train_ts <- ts(dat_train$value, frequency = 24)

aic_vec <- numeric(4)
bic_vec <- numeric(4)
adj_r2  <- numeric(4)

for (k in 1:4) {
  Xk <- as.data.frame(fourier(y_train_ts, K = k))
  colnames(Xk) <- make.names(colnames(Xk))
  mk <- lm(dat_train$value ~ ., data = Xk)
  aic_vec[k] <- AIC(mk)
  bic_vec[k] <- BIC(mk)
  adj_r2[k]  <- summary(mk)$adj.r.squared
  cat(sprintf("  K=%d:  AIC=%.2f  BIC=%.2f  Adj.R²=%.4f\n",
              k, aic_vec[k], bic_vec[k], adj_r2[k]))
}

best_k <- which.min(aic_vec)
cat(sprintf("Best K = %d  (AIC = %.2f)\n", best_k, aic_vec[best_k]))

t2 <- proc.time()
X_tr <- as.data.frame(fourier(y_train_ts, K = best_k))
colnames(X_tr) <- make.names(colnames(X_tr))
m2 <- lm(dat_train$value ~ ., data = X_tr)
t2_elapsed <- (proc.time() - t2)["elapsed"]

X_te <- as.data.frame(fourier(y_train_ts, K = best_k, h = n_test))
colnames(X_te) <- make.names(colnames(X_te))
pred_m2_train <- fitted(m2)
pred_m2_test  <- predict(m2, newdata = X_te)
resid_m2      <- residuals(m2)

m2_tr <- compute_metrics(dat_train$value, pred_m2_train)
m2_te <- compute_metrics(dat_test$value,  pred_m2_test)
lb_m2 <- Box.test(resid_m2, lag = 48, type = "Ljung-Box")$p.value

cat(sprintf("AIC (K=%d) = %.2f\n", best_k, AIC(m2)))
cat(sprintf("Time       = %.3f s\n", t2_elapsed))
cat(sprintf("Train:  RMSE=%.2f  MAE=%.2f  MAPE=%.2f%%\n",
            m2_tr$RMSE, m2_tr$MAE, m2_tr$MAPE))
cat(sprintf("Test:   RMSE=%.2f  MAE=%.2f  MAPE=%.2f%%\n",
            m2_te$RMSE, m2_te$MAE, m2_te$MAPE))
cat(sprintf("Ljung-Box Q(48) p = %.4f\n", lb_m2))

img("fig_s4_04_model2_k_selection.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  plot(1:4, aic_vec, type = "b", pch = 16, col = "steelblue", lwd = 2,
       main = "Mô hình 2 (Fourier): AIC theo số sóng K",
       xlab = "K (số sóng hài)", ylab = "AIC",
       xaxt = "n", cex.main = 1.0)
  axis(1, at = 1:4)
  points(best_k, aic_vec[best_k], col = "red", pch = 8, cex = 2.5, lwd = 2)
  legend("topright",
         c("AIC theo K", sprintf("K tốt nhất = %d", best_k)),
         col = c("steelblue", "red"), lty = c(1, NA),
         pch = c(16, 8), lwd = c(2, 2), cex = 0.9, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 700, height = 500)

img("fig_s4_05_model2_fit.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  plot(dat_test$datetime, dat_test$value, type = "l", col = "gray50",
       main = sprintf("Mô hình 2 (Fourier K=%d): Dự báo vs Thực tế", best_k),
       xlab = "Thời gian", ylab = "Tiêu thụ (MW)", cex.main = 1.0)
  lines(dat_test$datetime, pred_m2_test, col = "darkgreen", lwd = 1.5)
  legend("topright", c("Thực tế", "Dự báo M2"),
         col = c("gray50", "darkgreen"), lty = 1, lwd = c(1, 1.5),
         cex = 0.85, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 1200, height = 450)

img("fig_s4_06_model2_residuals.png", expression({
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  plot(dat_train$datetime, resid_m2, type = "l", col = "darkgreen",
       main = "Phần dư M2 theo thời gian",
       xlab = "Thời gian", ylab = "Phần dư (MW)", cex.main = 0.95)
  abline(h = 0, col = "red", lty = 2)
  acf(resid_m2, lag.max = 72, main = "ACF phần dư M2", cex.main = 0.95)
  qqnorm(resid_m2, main = "Q-Q phần dư M2", cex.main = 0.95)
  qqline(resid_m2, col = "red")
  par(mfrow = c(1, 1))
}), width = 1400, height = 500)

# =============================================================================
# 5. MODEL 3 — SARIMA
# =============================================================================
cat("\n=== MODEL 3: SARIMA ===\n")

t3 <- proc.time()
m3 <- auto.arima(y_train_ts, seasonal = TRUE, D = 1,
                 max.p = 2, max.q = 2, max.P = 1, max.Q = 1,
                 stepwise = TRUE, approximation = TRUE,
                 trace = FALSE)
t3_elapsed <- (proc.time() - t3)["elapsed"]

p3 <- m3$arma[1]; d3 <- m3$arma[6]; q3 <- m3$arma[2]
P3 <- m3$arma[3]; D3 <- m3$arma[7]; Q3 <- m3$arma[4]
cat(sprintf("SARIMA (%d,%d,%d)(%d,%d,%d)[24]\n", p3, d3, q3, P3, D3, Q3))
cat(sprintf("AIC  = %.2f\n", m3$aic))
cat(sprintf("Time = %.3f s\n", t3_elapsed))

fc3           <- forecast(m3, h = n_test)
pred_m3_test  <- as.numeric(fc3$mean)
pred_m3_train <- as.numeric(fitted(m3))
resid_m3      <- as.numeric(residuals(m3))

m3_tr <- compute_metrics(as.numeric(y_train_ts), pred_m3_train)
m3_te <- compute_metrics(dat_test$value, pred_m3_test)
lb_m3 <- Box.test(resid_m3, lag = 48, type = "Ljung-Box")$p.value

cat(sprintf("Train:  RMSE=%.2f  MAE=%.2f  MAPE=%.2f%%\n",
            m3_tr$RMSE, m3_tr$MAE, m3_tr$MAPE))
cat(sprintf("Test:   RMSE=%.2f  MAE=%.2f  MAPE=%.2f%%\n",
            m3_te$RMSE, m3_te$MAE, m3_te$MAPE))
cat(sprintf("Ljung-Box Q(48) p = %.4f\n", lb_m3))

img("fig_s4_07_model3_fit.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  plot(dat_test$datetime, dat_test$value, type = "l", col = "gray50",
       main = sprintf("Mô hình 3 SARIMA(%d,%d,%d)(%d,%d,%d)[24]: Dự báo vs Thực tế",
                      p3, d3, q3, P3, D3, Q3),
       xlab = "Thời gian", ylab = "Tiêu thụ (MW)", cex.main = 0.95)
  lines(dat_test$datetime, pred_m3_test, col = "darkorchid", lwd = 1.5)
  legend("topright", c("Thực tế", "Dự báo M3"),
         col = c("gray50", "darkorchid"), lty = 1, lwd = c(1, 1.5),
         cex = 0.85, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 1200, height = 450)

img("fig_s4_08_model3_residuals.png", expression({
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  plot(resid_m3, type = "l", col = "darkorchid",
       main = "Phần dư M3 theo thời gian",
       xlab = "Quan sát", ylab = "Phần dư (MW)", cex.main = 0.95)
  abline(h = 0, col = "red", lty = 2)
  acf(resid_m3, lag.max = 72, main = "ACF phần dư M3", cex.main = 0.95)
  qqnorm(resid_m3, main = "Q-Q phần dư M3", cex.main = 0.95)
  qqline(resid_m3, col = "red")
  par(mfrow = c(1, 1))
}), width = 1400, height = 500)

# =============================================================================
# 6. SO SÁNH TỔNG HỢP
# =============================================================================
cat("\n=== COMPARISON TABLE ===\n")
cat(sprintf("%-16s %10s %10s %10s %10s %10s\n",
            "Model", "AIC", "RMSE_te", "MAE_te", "MAPE_te%", "Time(s)"))
cat(sprintf("%-16s %10.2f %10.2f %10.2f %10.2f %10.3f\n",
            "M1 Indicator", AIC(m1), m1_te$RMSE, m1_te$MAE, m1_te$MAPE, t1_elapsed))
cat(sprintf("%-16s %10.2f %10.2f %10.2f %10.2f %10.3f\n",
            "M2 Fourier",   AIC(m2), m2_te$RMSE, m2_te$MAE, m2_te$MAPE, t2_elapsed))
cat(sprintf("%-16s %10.2f %10.2f %10.2f %10.2f %10.3f\n",
            "M3 SARIMA",   m3$aic,  m3_te$RMSE, m3_te$MAE, m3_te$MAPE, t3_elapsed))

# Bảng train metrics
cat("\n=== TRAIN METRICS ===\n")
cat(sprintf("%-16s %10s %10s %10s\n", "Model", "RMSE_tr", "MAE_tr", "MAPE_tr%"))
cat(sprintf("%-16s %10.2f %10.2f %10.2f\n", "M1 Indicator", m1_tr$RMSE, m1_tr$MAE, m1_tr$MAPE))
cat(sprintf("%-16s %10.2f %10.2f %10.2f\n", "M2 Fourier",   m2_tr$RMSE, m2_tr$MAE, m2_tr$MAPE))
cat(sprintf("%-16s %10.2f %10.2f %10.2f\n", "M3 SARIMA",    m3_tr$RMSE, m3_tr$MAE, m3_tr$MAPE))

# Ljung-Box summary
cat("\n=== LJUNG-BOX Q(48) p-values ===\n")
cat(sprintf("M1: %.4f\nM2: %.4f\nM3: %.4f\n", lb_m1, lb_m2, lb_m3))

# Fig so sánh toàn test
img("fig_s5_01_forecast_compare.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  plot(dat_test$datetime, dat_test$value, type = "l",
       col = "black", lwd = 0.8,
       main = "So sánh dự báo ba mô hình — tập kiểm tra",
       xlab = "Thời gian", ylab = "Tiêu thụ (MW)", cex.main = 1.0)
  lines(dat_test$datetime, pred_m1_test, col = "steelblue",  lwd = 1.5, lty = 1)
  lines(dat_test$datetime, pred_m2_test, col = "darkgreen",  lwd = 1.5, lty = 2)
  lines(dat_test$datetime, pred_m3_test, col = "darkorchid", lwd = 1.5, lty = 3)
  legend("topright",
         c("Thực tế", "M1 Indicator", "M2 Fourier", "M3 SARIMA"),
         col  = c("black", "steelblue", "darkgreen", "darkorchid"),
         lty  = c(1, 1, 2, 3),
         lwd  = c(0.8, 1.5, 1.5, 1.5),
         cex  = 0.85, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 1400, height = 500)

# Fig zoom 7 ngày đầu test
n_zoom <- min(168, n_test)
img("fig_s5_02_forecast_zoom.png", expression({
  par(mar = c(4, 4.5, 3, 1))
  plot(dat_test$datetime[1:n_zoom], dat_test$value[1:n_zoom],
       type = "l", col = "black", lwd = 1,
       main = "So sánh dự báo — 7 ngày đầu tập kiểm tra",
       xlab = "Thời gian", ylab = "Tiêu thụ (MW)", cex.main = 1.0)
  lines(dat_test$datetime[1:n_zoom], pred_m1_test[1:n_zoom],
        col = "steelblue",  lwd = 1.5, lty = 1)
  lines(dat_test$datetime[1:n_zoom], pred_m2_test[1:n_zoom],
        col = "darkgreen",  lwd = 1.5, lty = 2)
  lines(dat_test$datetime[1:n_zoom], pred_m3_test[1:n_zoom],
        col = "darkorchid", lwd = 1.5, lty = 3)
  abline(v = dat_test$datetime[seq(1, n_zoom, 24)],
         col = "gray70", lty = 2, lwd = 0.5)
  legend("topright",
         c("Thực tế", "M1 Indicator", "M2 Fourier", "M3 SARIMA"),
         col  = c("black", "steelblue", "darkgreen", "darkorchid"),
         lty  = c(1, 1, 2, 3),
         lwd  = c(1, 1.5, 1.5, 1.5),
         cex  = 0.85, bg = "white")
  grid(col = "gray90", lty = 1)
}), width = 1200, height = 450)

# Fig ACF phần dư so sánh
img("fig_s5_03_residual_acf_compare.png", expression({
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  acf(resid_m1, lag.max = 72,
      main = "ACF phần dư M1 (Indicator)", cex.main = 0.9)
  acf(resid_m2, lag.max = 72,
      main = "ACF phần dư M2 (Fourier)",   cex.main = 0.9)
  acf(resid_m3, lag.max = 72,
      main = "ACF phần dư M3 (SARIMA)",    cex.main = 0.9)
  par(mfrow = c(1, 1))
}), width = 1400, height = 500)

cat("\n✓ Hoàn tất models_report.R — ảnh đã lưu vào report/img/\n")
