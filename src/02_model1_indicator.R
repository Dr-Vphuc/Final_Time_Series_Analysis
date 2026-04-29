# =============================================================================
# 02_model1_indicator.R — Mô hình 1: Hồi quy mùa vụ chỉ số (24 hệ số giờ)
# Tương ứng NGÀY 2 trong timeline.
# =============================================================================

source("src/00_setup.R")

# --- 1. Load dữ liệu đã xử lý ------------------------------------------------
processed <- readRDS(file.path(DATA_DIR, "processed_data.rds"))
y   <- processed$y
dat <- processed$dat

# --- 2. Chia train/test ------------------------------------------------------
split_obj <- train_test_split(as.numeric(y), train_ratio = 0.8)
y_train <- split_obj$train
y_test  <- split_obj$test
n_train <- split_obj$n_train
n_test  <- split_obj$n_test

train_dat <- dat[1:n_train, ]
test_dat  <- dat[(n_train + 1):nrow(dat), ]

cat("Train:", n_train, "obs | Test:", n_test, "obs\n\n")

# --- 3. Tạo biến cho mô hình -------------------------------------------------
train_dat$t          <- 1:n_train
train_dat$hour_fact  <- factor(train_dat$hour, levels = 0:23)
train_dat$is_weekend <- as.integer(wday(train_dat$datetime) %in% c(1, 7))

test_dat$t          <- (n_train + 1):(n_train + n_test)
test_dat$hour_fact  <- factor(test_dat$hour, levels = 0:23)
test_dat$is_weekend <- as.integer(wday(test_dat$datetime) %in% c(1, 7))

# --- 4. Fit Model 1 ----------------------------------------------------------
cat("=== FIT MODEL 1: Hồi quy mùa chỉ số ===\n")
model1 <- lm(value ~ hour_fact + is_weekend + t, data = train_dat)
print(summary(model1))

cat("\nAIC(model1) =", AIC(model1), "\n")
cat("BIC(model1) =", BIC(model1), "\n\n")

# --- 5. Chẩn đoán phần dư ----------------------------------------------------
res1 <- rstudent(model1)

save_plot("02a_model1_residuals.png", expression({
  par(mfrow = c(2, 2))
  plot(res1, type = "l", main = "M1: Phần dư chuẩn hóa", ylab = "rstudent")
  abline(h = c(-2, 0, 2), lty = 2, col = "red")
  qqnorm(res1, main = "M1: QQ-plot")
  qqline(res1, col = "red")
  acf(res1, lag.max = 72, main = "M1: ACF phần dư (kỳ vọng: còn nhiều)")
  hist(res1, breaks = 50, main = "M1: Histogram phần dư", col = "lightblue")
  par(mfrow = c(1, 1))
}))

cat("[Shapiro-Wilk] (mẫu nhỏ — chỉ chạy nếu n < 5000)\n")
if (length(res1) <= 5000) {
  print(shapiro.test(res1))
} else {
  cat("  → mẫu quá lớn, dùng QQ-plot để đánh giá.\n")
}

cat("\n[Ljung-Box] tại lag = 24, 48, 72\n")
for (lg in c(24, 48, 72)) {
  bt <- Box.test(res1, lag = lg, type = "Ljung-Box")
  cat(sprintf("  lag = %2d | p-value = %.4g\n", lg, bt$p.value))
}

# --- 6. Dự báo trên test -----------------------------------------------------
pred1 <- predict(model1, newdata = test_dat)

metrics1 <- compute_metrics(y_test, pred1)
metrics1$Model <- "M1_Indicator"
metrics1$AIC   <- AIC(model1)
cat("\n=== METRICS MODEL 1 ===\n"); print(metrics1)

# --- 7. Vẽ dự báo ------------------------------------------------------------
save_plot("02b_model1_forecast.png", expression({
  plot(y_test, type = "l", col = "black", lwd = 1.5,
       main = "Model 1: Actual vs Predicted (test)",
       ylab = "MW", xlab = "Giờ trong test")
  lines(pred1, col = "red", lwd = 1.5)
  legend("topright", c("Actual", "Predicted M1"),
         col = c("black", "red"), lty = 1, lwd = 1.5)
}))

# --- 8. Lưu kết quả ----------------------------------------------------------
saveRDS(list(model = model1, pred = pred1, metrics = metrics1, res = res1),
        file = file.path(OUTPUT_DIR, "model1_result.rds"))

cat("\n✓ Hoàn tất Model 1. Kết quả lưu trong output/model1_result.rds\n")
cat("  Tiếp theo: chạy 03_model2_fourier.R\n")
