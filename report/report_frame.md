1. Giới thiệu
   1.1. Vì sao nghiên cứu chuỗi tiêu thụ điện
   1.2. Mục tiêu: so sánh 3 phương pháp xử lý mùa vụ

2. Dữ liệu & EDA
   2.1. Nguồn dữ liệu, mô tả
   2.2. Trực quan hóa: chu kỳ ngày, chu kỳ tuần, xu hướng
   2.3. Kiểm định ADF / KPSS

3. Cơ sở lý thuyết
   3.1. Hồi quy mùa vụ chỉ số
   3.2. Hồi quy điều hòa (Fourier)
   3.3. SARIMA
   3.4. Tiêu chí so sánh: AIC, RMSE, MAE, Ljung-Box

4. Triển khai
   4.1. Tiền xử lý
   4.2. Chia tập huấn luyện / kiểm tra
   4.3. Mô hình 1: hồi quy mùa chỉ số
   4.4. Mô hình 2: hồi quy điều hòa (K = 1..4)
   4.5. Mô hình 3: SARIMA(p,d,q)(P,D,Q)_24

5. Kết quả & So sánh
   5.1. Bảng tổng hợp AIC, RMSE, MAE
   5.2. So sánh trực quan dự báo
   5.3. So sánh phần dư
   5.4. Thời gian tính toán

6. Thảo luận
   6.1. Khi nào nên dùng phương pháp nào (đối chiếu lý thuyết)
   6.2. Hạn chế của project
   6.3. Hướng mở rộng (nhiều chu kỳ, TBATS)

7. Kết luận
