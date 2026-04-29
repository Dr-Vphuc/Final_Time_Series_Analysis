# Project: So sánh 3 phương pháp xử lý mùa vụ trên dữ liệu tiêu thụ điện

Bộ template R code khung sẵn cho project cuối kỳ môn Phân tích chuỗi thời gian.

## Cấu trúc

```
src/
├── 00_setup.R              # Cài package, hàm tiện ích, đường dẫn
├── 01_data_eda.R           # Ngày 1: Tải dữ liệu + EDA + kiểm định tính dừng
├── 02_model1_indicator.R   # Ngày 2: Mô hình 1 — Hồi quy mùa chỉ số
├── 03_model2_fourier.R     # Ngày 3: Mô hình 2 — Hồi quy điều hòa Fourier (K=1..6)
├── 04_model3_sarima.R      # Ngày 4: Mô hình 3 — SARIMA
└── 05_compare.R            # Ngày 5: So sánh 3 mô hình + tổng hợp
```

Output sẽ tự động ghi vào:
- `data/`   — dữ liệu thô + file `.rds` đã xử lý
- `output/` — biểu đồ PNG, bảng CSV, kết quả model `.rds`

## Chuẩn bị

1. **Tải dataset** từ Kaggle: [PJM Hourly Energy Consumption](https://www.kaggle.com/datasets/robikscube/hourly-energy-consumption)
2. Đặt file `AEP_hourly.csv` vào thư mục `data/` ở thư mục gốc project.
3. Mở R / RStudio, đặt working directory tại thư mục gốc:
   ```r
   setwd("f:/Learning/Time_Series_Analysis")
   ```

## Cách chạy

Chạy lần lượt từng file (kết quả mỗi file được lưu lại để file sau dùng):

```r
source("src/01_data_eda.R")           # Ngày 1
source("src/02_model1_indicator.R")   # Ngày 2
source("src/03_model2_fourier.R")     # Ngày 3
source("src/04_model3_sarima.R")      # Ngày 4 — chậm nhất, ~10-30 phút
source("src/05_compare.R")            # Ngày 5
```

## Tinh chỉnh nhanh

- **Đổi phạm vi dữ liệu**: sửa `START_DATE` và `END_DATE` trong [01_data_eda.R](01_data_eda.R)
- **Đổi tỷ lệ train/test**: sửa `train_ratio` trong các file 02-04 (mặc định 0.8)
- **Đổi K Fourier**: sửa vòng lặp `for (K in 1:6)` trong [03_model2_fourier.R](03_model2_fourier.R)
- **Đổi ứng viên SARIMA**: sửa các dòng `fit_sarima_safe(...)` trong [04_model3_sarima.R](04_model3_sarima.R)
- **Bỏ qua SARIMA tay** (chỉ dùng `auto.arima`): comment ứng viên c1, c2 trong file 04

## Khi gặp lỗi

| Lỗi | Cách xử lý |
|-----|-----------|
| Không tìm thấy `AEP_hourly.csv` | Tải về và đặt vào `data/` |
| SARIMA fit quá lâu | Giảm phạm vi dữ liệu (1-2 tháng), hoặc chỉ dùng `auto.arima` |
| `auto.arima` không hội tụ | Đặt `approximation = TRUE`, `stepwise = TRUE` (mặc định) |
| Out of memory | Giảm số quan sát; đổi `frequency = 24` thay vì `168` |

## Kết quả mong đợi

Sau khi chạy xong toàn bộ:

- `output/compare_table.csv` — bảng AIC / RMSE / MAE / MAPE của 3 mô hình
- `output/ljungbox_table.csv` — p-value Ljung-Box ở các lag
- 12+ biểu đồ PNG để chèn vào báo cáo & slide
- 3 file `model{1,2,3}_result.rds` chứa kết quả đầy đủ (model, dự báo, phần dư)

Dùng các file này làm tư liệu trực tiếp cho báo cáo + slide trình bày.
