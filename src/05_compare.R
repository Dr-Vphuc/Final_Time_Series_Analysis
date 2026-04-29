# =============================================================================
# 05_compare.R — Tổng hợp so sánh 3 mô hình
# Tương ứng NGÀY 5 trong timeline.
# =============================================================================

source("src/00_setup.R")

# --- 1. Load kết quả 3 mô hình -----------------------------------------------
m1 <- readRDS(file.path(OUTPUT_DIR, "model1_result.rds"))
m2 <- readRDS(file.path(OUTPUT_DIR, "model2_result.rds"))
m3 <- readRDS(file.path(OUTPUT_DIR, "model3_result.rds"))

processed <- readRDS(file.path(DATA_DIR, "processed_data.rds"))
y         <- processed$y
split_obj <- train_test_split(as.numeric(y), train_ratio = 0.8)
y_test    <- split_obj$test
n_test    <- split_obj$n_test

# --- 2. Bảng so sánh tổng ----------------------------------------------------
compare_table <- rbind(m1$metrics, m2$metrics, m3$metrics)
compare_table <- compare_table[, c("Model", "AIC", "RMSE", "MAE", "MAPE")]
cat("\n========== BẢNG SO SÁNH 3 MÔ HÌNH ==========\n")
print(compare_table, row.names = FALSE)

# Lưu CSV
write.csv(compare_table,
          file = file.path(OUTPUT_DIR, "compare_table.csv"),
          row.names = FALSE)

# --- 3. Bảng Ljung-Box --------------------------------------------------------
lb_test <- function(res, lag, fitdf = 0) {
  bt <- Box.test(res, lag = lag, type = "Ljung-Box", fitdf = fitdf)
  bt$p.value
}

lb_table <- data.frame(
  Model      = c("M1_Indicator", paste0("M2_Fourier_K", m2$K_opt),
                 "M3_SARIMA"),
  LB_lag24   = c(lb_test(m1$res, 24),
                 lb_test(m2$res, 24),
                 lb_test(m3$res, 24, fitdf = length(coef(m3$model)))),
  LB_lag48   = c(lb_test(m1$res, 48),
                 lb_test(m2$res, 48),
                 lb_test(m3$res, 48, fitdf = length(coef(m3$model)))),
  LB_lag72   = c(lb_test(m1$res, 72),
                 lb_test(m2$res, 72),
                 lb_test(m3$res, 72, fitdf = length(coef(m3$model))))
)
cat("\n========== p-VALUE LJUNG-BOX ==========\n")
cat("(p > 0.05 => phần dư là nhiễu trắng)\n\n")
print(lb_table, row.names = FALSE)

write.csv(lb_table,
          file = file.path(OUTPUT_DIR, "ljungbox_table.csv"),
          row.names = FALSE)

# --- 4. Biểu đồ so sánh tổng — actual vs 3 mô hình ---------------------------
save_plot("05a_compare_forecasts.png", expression({
  plot(y_test, type = "l", col = "black", lwd = 2,
       main = "So sánh dự báo 3 mô hình trên tập test",
       ylab = "MW", xlab = "Giờ trong test")
  lines(m1$pred, col = "red",       lwd = 1.2)
  lines(m2$pred, col = "blue",      lwd = 1.2)
  lines(m3$pred, col = "darkgreen", lwd = 1.2)
  legend("topright",
         c("Actual", "M1 Indicator", paste0("M2 Fourier K=", m2$K_opt),
           "M3 SARIMA"),
         col = c("black", "red", "blue", "darkgreen"),
         lty = 1, lwd = c(2, 1.2, 1.2, 1.2))
}), width = 1200, height = 600)

# Zoom 1 tuần đầu của test
save_plot("05b_compare_forecasts_zoom.png", expression({
  zoom <- 1:min(168, n_test)
  plot(y_test[zoom], type = "l", col = "black", lwd = 2,
       main = "So sánh dự báo (zoom 1 tuần đầu test)",
       ylab = "MW", xlab = "Giờ")
  lines(zoom, m1$pred[zoom], col = "red",       lwd = 1.2)
  lines(zoom, m2$pred[zoom], col = "blue",      lwd = 1.2)
  lines(zoom, m3$pred[zoom], col = "darkgreen", lwd = 1.2)
  legend("topright",
         c("Actual", "M1", "M2", "M3"),
         col = c("black", "red", "blue", "darkgreen"),
         lty = 1, lwd = c(2, 1.2, 1.2, 1.2))
}), width = 1200, height = 600)

# --- 5. Biểu đồ ACF phần dư của 3 mô hình -----------------------------------
save_plot("05c_compare_residual_acf.png", expression({
  par(mfrow = c(3, 1))
  acf(m1$res, lag.max = 72, main = "ACF phần dư M1 (Indicator)")
  acf(m2$res, lag.max = 72, main = paste0("ACF phần dư M2 (Fourier K=", m2$K_opt, ")"))
  acf(m3$res, lag.max = 72, main = "ACF phần dư M3 (SARIMA)")
  par(mfrow = c(1, 1))
}), height = 900)

# --- 6. Bảng RMSE chi tiết theo giờ trong ngày ------------------------------
hour_test <- 0:(n_test - 1) %% 24
rmse_by_hour <- data.frame(
  hour = 0:23,
  M1   = sapply(0:23, function(h) sqrt(mean((y_test[hour_test == h] - m1$pred[hour_test == h])^2))),
  M2   = sapply(0:23, function(h) sqrt(mean((y_test[hour_test == h] - m2$pred[hour_test == h])^2))),
  M3   = sapply(0:23, function(h) sqrt(mean((y_test[hour_test == h] - m3$pred[hour_test == h])^2)))
)

save_plot("05d_rmse_by_hour.png", expression({
  matplot(rmse_by_hour$hour,
          rmse_by_hour[, c("M1", "M2", "M3")],
          type = "b", lty = 1, pch = 16,
          col = c("red", "blue", "darkgreen"),
          xlab = "Giờ trong ngày", ylab = "RMSE",
          main = "RMSE phân theo giờ trong ngày")
  legend("topleft", c("M1", "M2", "M3"),
         col = c("red", "blue", "darkgreen"), lty = 1, pch = 16)
}))

# --- 7. In kết luận sơ bộ ---------------------------------------------------
best_rmse <- compare_table$Model[which.min(compare_table$RMSE)]
best_aic  <- compare_table$Model[which.min(compare_table$AIC)]

cat("\n========== KẾT LUẬN SƠ BỘ ==========\n")
cat("• Mô hình tốt nhất theo RMSE :", best_rmse, "\n")
cat("• Mô hình tốt nhất theo AIC  :", best_aic,  "\n")
cat("\nXem các file PNG trong output/ để dùng cho báo cáo & slide.\n")
cat("File CSV: compare_table.csv, ljungbox_table.csv\n")
