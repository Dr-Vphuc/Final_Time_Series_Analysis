# =============================================================================
# 03_model2_fourier.R — Mô hình 2: Hồi quy điều hòa Fourier (K = 1..6)
# Tương ứng NGÀY 3 trong timeline.
# =============================================================================

source("src/00_setup.R")

# --- 1. Load dữ liệu ---------------------------------------------------------
processed <- readRDS(file.path(DATA_DIR, "processed_data.rds"))
y   <- processed$y
dat <- processed$dat

split_obj <- train_test_split(as.numeric(y), train_ratio = 0.8)
y_train   <- split_obj$train
y_test    <- split_obj$test
n_train   <- split_obj$n_train
n_test    <- split_obj$n_test

train_dat <- dat[1:n_train, ]
test_dat  <- dat[(n_train + 1):nrow(dat), ]

train_dat$t          <- 1:n_train
test_dat$t           <- (n_train + 1):(n_train + n_test)
train_dat$is_weekend <- as.integer(wday(train_dat$datetime) %in% c(1, 7))
test_dat$is_weekend  <- as.integer(wday(test_dat$datetime)  %in% c(1, 7))

# Tạo ts cho hàm fourier
y_train_ts <- ts(y_train, frequency = 24)

# --- 2. Hàm fit Fourier với K cho trước --------------------------------------
fit_fourier <- function(K) {
  X_train <- fourier(y_train_ts, K = K)
  X_test  <- fourier(y_train_ts, K = K, h = n_test)

  df_train <- cbind(
    data.frame(value = y_train, t = train_dat$t, is_weekend = train_dat$is_weekend),
    as.data.frame(X_train)
  )
  df_test <- cbind(
    data.frame(t = test_dat$t, is_weekend = test_dat$is_weekend),
    as.data.frame(X_test)
  )

  # Tên cột Fourier có dấu / hyphen — sửa để dùng được trong formula
  fnames <- colnames(X_train)
  fnames_clean <- gsub("[^[:alnum:]_]", "_", fnames)
  colnames(df_train)[(ncol(df_train) - length(fnames) + 1):ncol(df_train)] <- fnames_clean
  colnames(df_test)[(ncol(df_test)  - length(fnames) + 1):ncol(df_test)]   <- fnames_clean

  fmla <- as.formula(
    paste("value ~ t + is_weekend +", paste(fnames_clean, collapse = " + "))
  )
  model <- lm(fmla, data = df_train)

  pred  <- predict(model, newdata = df_test)
  list(K = K, model = model, pred = pred,
       AIC = AIC(model), BIC = BIC(model))
}

# --- 3. Quét K = 1..6 --------------------------------------------------------
cat("=== QUÉT K = 1..6 ===\n")
results <- list()
summary_table <- data.frame()

for (K in 1:6) {
  cat("Đang fit K =", K, "...\n")
  res <- fit_fourier(K)
  results[[K]] <- res

  m <- compute_metrics(y_test, res$pred)
  summary_table <- rbind(summary_table, data.frame(
    K = K, AIC = res$AIC, BIC = res$BIC,
    RMSE = m$RMSE, MAE = m$MAE, MAPE = m$MAPE
  ))
}

cat("\n=== BẢNG SO SÁNH K ===\n")
print(summary_table)

# --- 4. Chọn K tối ưu --------------------------------------------------------
best_K_aic  <- summary_table$K[which.min(summary_table$AIC)]
best_K_rmse <- summary_table$K[which.min(summary_table$RMSE)]
cat("\n→ K tối ưu theo AIC :", best_K_aic, "\n")
cat("→ K tối ưu theo RMSE:", best_K_rmse, "\n")

K_OPT <- best_K_aic
best  <- results[[K_OPT]]

# --- 5. Chẩn đoán phần dư cho K tối ưu ---------------------------------------
res2 <- rstudent(best$model)

save_plot("03a_model2_residuals.png", expression({
  par(mfrow = c(2, 2))
  plot(res2, type = "l", main = paste0("M2 (K=", K_OPT, "): Phần dư chuẩn hóa"))
  abline(h = c(-2, 0, 2), lty = 2, col = "red")
  qqnorm(res2, main = paste0("M2 (K=", K_OPT, "): QQ-plot"))
  qqline(res2, col = "red")
  acf(res2, lag.max = 72, main = "M2: ACF phần dư")
  hist(res2, breaks = 50, main = "M2: Histogram phần dư", col = "lightyellow")
  par(mfrow = c(1, 1))
}))

cat("\n[Ljung-Box] cho M2 tại lag = 24, 48, 72\n")
for (lg in c(24, 48, 72)) {
  bt <- Box.test(res2, lag = lg, type = "Ljung-Box")
  cat(sprintf("  lag = %2d | p-value = %.4g\n", lg, bt$p.value))
}

# --- 6. Vẽ so sánh thành phần mùa M1 vs M2 -----------------------------------
m1_res <- readRDS(file.path(OUTPUT_DIR, "model1_result.rds"))
m1_coef <- coef(m1_res$model)
m1_hour_effect <- c(0, m1_coef[grepl("hour_fact", names(m1_coef))])

# Thành phần mùa của M2: chỉ giữ phần Fourier, đặt t = 0, weekend = 0
X_one_day <- fourier(y_train_ts, K = K_OPT)[1:24, , drop = FALSE]
m2_coef <- coef(best$model)
fourier_idx <- grep("S|C", names(m2_coef))   # các cột Fourier trong coef
m2_hour_effect <- as.numeric(X_one_day %*% m2_coef[fourier_idx])

save_plot("03b_compare_seasonal_M1_M2.png", expression({
  plot(0:23, m1_hour_effect, type = "s", col = "red", lwd = 2,
       main = "Thành phần mùa: M1 (bậc thang) vs M2 (mượt)",
       xlab = "Giờ", ylab = "Hiệu ứng tương đối (MW)",
       ylim = range(c(m1_hour_effect, m2_hour_effect)))
  lines(0:23, m2_hour_effect, col = "blue", lwd = 2)
  legend("topright", c("M1: Indicator", paste0("M2: Fourier K=", K_OPT)),
         col = c("red", "blue"), lty = 1, lwd = 2)
}))

# --- 7. Vẽ dự báo M2 trên test -----------------------------------------------
save_plot("03c_model2_forecast.png", expression({
  plot(y_test, type = "l", col = "black", lwd = 1.5,
       main = paste0("Model 2 (K=", K_OPT, "): Actual vs Predicted (test)"),
       ylab = "MW", xlab = "Giờ trong test")
  lines(best$pred, col = "blue", lwd = 1.5)
  legend("topright", c("Actual", paste0("Predicted M2 (K=", K_OPT, ")")),
         col = c("black", "blue"), lty = 1, lwd = 1.5)
}))

# --- 8. Lưu kết quả ----------------------------------------------------------
metrics2 <- compute_metrics(y_test, best$pred)
metrics2$Model <- paste0("M2_Fourier_K", K_OPT)
metrics2$AIC   <- best$AIC

saveRDS(list(K_opt = K_OPT, model = best$model, pred = best$pred,
             metrics = metrics2, res = res2, scan_table = summary_table),
        file = file.path(OUTPUT_DIR, "model2_result.rds"))

cat("\n✓ Hoàn tất Model 2. K tối ưu =", K_OPT, "\n")
cat("  Tiếp theo: chạy 04_model3_sarima.R\n")
