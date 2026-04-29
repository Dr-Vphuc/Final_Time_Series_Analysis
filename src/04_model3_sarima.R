# =============================================================================
# 04_model3_sarima.R — Mô hình 3: SARIMA(p,d,q)(P,D,Q)_24
# Tương ứng NGÀY 4 trong timeline (ngày dài nhất).
#
# CẢNH BÁO: SARIMA fit có thể chậm 10-30 phút.
# Nếu không hội tụ, dùng auto.arima() làm fallback.
# =============================================================================

source("src/00_setup.R")

# --- 1. Load dữ liệu ---------------------------------------------------------
processed <- readRDS(file.path(DATA_DIR, "processed_data.rds"))
y   <- processed$y
dat <- processed$dat

split_obj <- train_test_split(as.numeric(y), train_ratio = 0.8)
y_train   <- ts(split_obj$train, frequency = 24)
y_test    <- split_obj$test
n_test    <- split_obj$n_test

cat("Train:", length(y_train), " | Test:", n_test, "\n\n")

# --- 2. Phân tích sơ bộ trước khi fit ----------------------------------------
cat("=== KIỂM TRA TÍNH DỪNG TRÊN Y_TRAIN ===\n")
cat("[ADF]\n"); print(adf.test(y_train))

cat("\n[ADF sau sai phân mùa lag=24]\n")
y_diff_seas <- diff(y_train, lag = 24)
print(adf.test(y_diff_seas))

save_plot("04a_acf_after_seasonal_diff.png", expression({
  par(mfrow = c(2, 1))
  acf(y_diff_seas,  lag.max = 72, main = "ACF sau sai phân mùa")
  pacf(y_diff_seas, lag.max = 72, main = "PACF sau sai phân mùa")
  par(mfrow = c(1, 1))
}))

# EACF cho phần không mùa (chỉ chạy nếu mẫu không quá lớn)
if (length(y_diff_seas) <= 5000) {
  cat("\n[EACF cho phần không mùa]\n")
  print(eacf(y_diff_seas, ar.max = 5, ma.max = 5))
}

# --- 3. Fit auto.arima làm baseline (nhanh) ----------------------------------
cat("\n=== FIT AUTO.ARIMA BASELINE ===\n")
cat("(có thể mất vài phút...)\n")
t0 <- Sys.time()
model_auto <- auto.arima(
  y_train,
  seasonal = TRUE,
  stepwise = TRUE,         # đặt FALSE để quét đầy đủ (chậm hơn nhiều)
  approximation = TRUE,    # đặt FALSE để chính xác hơn (chậm hơn)
  trace = FALSE
)
cat("Thời gian:", round(difftime(Sys.time(), t0, units = "secs"), 1), "giây\n")
cat("\nMô hình auto.arima chọn:\n"); print(model_auto)

# --- 4. Fit thủ công 2 ứng viên SARIMA (tùy chọn — bỏ qua nếu chậm) ----------
fit_sarima_safe <- function(p, d, q, P, D, Q, S = 24, label = "") {
  cat("\n[Đang fit SARIMA(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")_", S,
      "] ", label, "\n", sep = "")
  t0 <- Sys.time()
  fit <- tryCatch(
    Arima(y_train, order = c(p, d, q),
          seasonal = list(order = c(P, D, Q), period = S),
          method = "CSS-ML"),
    error = function(e) { cat("LỖI:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(fit)) {
    cat("  → AIC =", AIC(fit), "| time =",
        round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
  }
  fit
}

candidates <- list()

# Ứng viên 1: dựa trên ACF/PACF gợi ý đơn giản
candidates[["c1_(1,0,1)(1,1,1)"]] <- fit_sarima_safe(1, 0, 1, 1, 1, 1)

# Ứng viên 2: ít tham số hơn để so sánh
candidates[["c2_(2,0,0)(0,1,1)"]] <- fit_sarima_safe(2, 0, 0, 0, 1, 1)

# Ứng viên 3: từ auto.arima
candidates[["c3_auto"]] <- model_auto

# --- 5. Chọn mô hình tốt nhất theo AIC ---------------------------------------
candidates <- candidates[!sapply(candidates, is.null)]
aic_table  <- data.frame(
  model = names(candidates),
  AIC   = sapply(candidates, AIC),
  BIC   = sapply(candidates, BIC)
)
cat("\n=== BẢNG AIC CÁC ỨNG VIÊN ===\n"); print(aic_table)

best_name  <- aic_table$model[which.min(aic_table$AIC)]
model3     <- candidates[[best_name]]
cat("\n→ Chọn:", best_name, "\n")
print(model3)

# --- 6. Chẩn đoán phần dư ----------------------------------------------------
res3 <- residuals(model3)

save_plot("04b_model3_residuals.png", expression({
  par(mfrow = c(2, 2))
  plot(res3, main = "M3: Phần dư", ylab = "")
  abline(h = 0, lty = 2, col = "red")
  qqnorm(res3, main = "M3: QQ-plot")
  qqline(res3, col = "red")
  acf(res3, lag.max = 72, main = "M3: ACF phần dư")
  hist(res3, breaks = 50, main = "M3: Histogram phần dư", col = "lightcoral")
  par(mfrow = c(1, 1))
}))

cat("\n[Ljung-Box] cho M3 tại lag = 24, 48, 72\n")
df_adj <- length(model3$coef)   # bậc tự do điều chỉnh
for (lg in c(24, 48, 72)) {
  bt <- Box.test(res3, lag = lg, type = "Ljung-Box", fitdf = df_adj)
  cat(sprintf("  lag = %2d | p-value = %.4g\n", lg, bt$p.value))
}

# --- 7. Dự báo trên test -----------------------------------------------------
cat("\n=== DỰ BÁO ===\n")
fc <- forecast(model3, h = n_test)
pred3 <- as.numeric(fc$mean)

metrics3 <- compute_metrics(y_test, pred3)
metrics3$Model <- paste0("M3_SARIMA_", best_name)
metrics3$AIC   <- AIC(model3)
cat("\nMetrics M3:\n"); print(metrics3)

# --- 8. Vẽ dự báo có khoảng tin cậy 95% --------------------------------------
save_plot("04c_model3_forecast.png", expression({
  plot(y_test, type = "l", col = "black", lwd = 1.5,
       main = paste0("Model 3 (", best_name, "): Actual vs Predicted (test)"),
       ylab = "MW", xlab = "Giờ trong test",
       ylim = range(c(y_test, fc$lower[, 2], fc$upper[, 2]), na.rm = TRUE))
  polygon(c(1:n_test, rev(1:n_test)),
          c(fc$lower[, 2], rev(fc$upper[, 2])),
          col = rgb(0, 1, 0, 0.2), border = NA)
  lines(pred3, col = "darkgreen", lwd = 1.5)
  legend("topright",
         c("Actual", "Predicted M3", "95% CI"),
         col = c("black", "darkgreen", rgb(0, 1, 0, 0.4)),
         lty = c(1, 1, NA), lwd = c(1.5, 1.5, NA),
         pch = c(NA, NA, 15))
}))

# --- 9. Lưu kết quả ----------------------------------------------------------
saveRDS(list(model = model3, pred = pred3, fc = fc, metrics = metrics3,
             res = res3, candidates_aic = aic_table),
        file = file.path(OUTPUT_DIR, "model3_result.rds"))

cat("\n✓ Hoàn tất Model 3. Kết quả lưu trong output/model3_result.rds\n")
cat("  Tiếp theo: chạy 05_compare.R\n")
