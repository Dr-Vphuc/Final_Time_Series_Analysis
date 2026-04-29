#18/3/2026

# QUÁ TRÌNH TRUNG BÌNH TRƯỢT TỰ HỒI QUY ARMA(p, q)

# x_t = phi_1*x_{t-1} + ... + phi_p*x_{t-p} + w_t + theta_1*w_{t-1} + ... + theta_q*w_{t-q}

# Đa thức AR
phi(z) = 1 - phi_1*z - ... - phi_p*z^p

# Đa thức MA
theta(z) = 1 + theta_1*z + ... + theta_q*z^q

# Chú ý:
# 1. Các quá trình ARMA(p, q) yêu cầu các đa thức AR và MA không có thừa số chung
# 2. Đối với các quá trình ARMA(p, q) tổng quát,
#    yêu cầu thỏa mãn cả tính dừng và tính khả nghịch

# Nhắc lại:
# + Tính dừng: Tất cả nghiệm của đa thức AR phải nằm ngoài đường tròn đơn vị
# + Tính khả nghịch: Tất cả nghiệm của đa thức MA phải nằm ngoài đường tròn đơn vị
#
#nên khó xác định p và q vì ACF và PACF giảm dần
# => dùng EACF để hỗ trợ

library(tseries)
library(astsa)

# Sinh dữ liệu theo mô hình ARMA(1,1) với tham số 
set.seed(123)
data <- sarima.sim(ar = 0.6, ma = 0.4, n = 1000)

# Vẽ chuỗi
plot(data)

# Kiểm định ADF (p-value < 0.05 => chuỗi dừng)
adf.test(data)

# Vẽ ACF và PACF
par(mfrow = c(1, 2))
acf(data)
pacf(data)

# Nhận thấy
# Các giá trị ACF giảm dần
# Các giá trị PACF giảm dần => gợi ý mô hình ARMA

# Ngoài ra, giá trị PACF tại 3 độ trễ đầu tiên khác 0 một cách đáng kể
# => có thể xét thêm các mô hình AR(2) hoặc AR(3)


# Để xác định bậc của mô hình ARMA thì ta cần xét bảng EACF
# Trong bảng này:
#   x: hệ số tương quan mở rộng có ý nghĩa thống kê
#   0: hệ số tương quan mở rộng không có ý nghĩa thống kê

# Ta cần xác định một vị trí (p, q) sao cho:
#   từ vị trí đó trở đi (bên phải và phía dưới) đều là 0
#   (tạo thành một tam giác số 0)
install.packages("TSA")
library(TSA)
eacf(data)

# Dựa vào biểu đồ EACF, ta dự đoán bậc của mô hình ARMA là ARMA(1,1)
# Ước lượng tham số và chẩn đoán mô hình
mod1 <- sarima(data, 2, 0, 0, no.constant = TRUE)   # mô hình AR(2)

# Các hệ số đều có p-value < 0.05 => đều có ý nghĩa thống kê

# Phương trình:
# x_t = 0.9361*x_{t-1} - 0.2791*x_{t-2} + w_t


# Chẩn đoán mô hình

# - Đồ thị phần dư chuẩn hóa: không có xu hướng
# - Không có dấu hiệu phương sai thay đổi
# - Đồ thị QQ-plot: các điểm dữ liệu rất sát đường thẳng chuẩn
# - Đồ thị ACF: không còn tự tương quan đáng kể

# đồ thị p_value của  Kiểm định Ljung-Box
# Kết quả:
# thấy rằng các giá trị p-value hầu hết đều < 0.05 ở các độ trễ

# => nghĩa là phần dư vẫn còn tự tương quan tại các độ trễ nhỏ
# => phần dư chưa hoàn toàn là nhiễu trắng

# Ghi chú: kiểm định Ljung-Box (tổng quát)
# H0: không còn tự tương quan trong dữ liệu

# Mô hình này tương đương với ARMA(1, 1)
mod3 <- sarima(data, 1, 0, 1, no.constant = TRUE)

# Dòng 85-86: Chú thích về hệ số và phương trình (dựa trên kết quả chạy)
# Các hệ số đều có p-value < 0.05 => đều có ý nghĩa thống kê
# Phương trình: x_t = 0.5609*x_{t-1} + w_t + 0.4076*w_{t-1}

# ---------------------------------------------------------
# Dòng 87: Chẩn đoán mô hình (Diagnostics)
# Khi chạy lệnh sarima(), R sẽ tự động xuất ra 4 đồ thị chẩn đoán:

# 1. Đồ thị phần dư chuẩn hóa (Standardized Residuals):
# => Không có xu hướng (no trend), không có dấu hiệu về phương sai thay đổi.

# 2. Đồ thị QQ-plot:
# => Các điểm dữ liệu rất sát đường thẳng, chứng tỏ phần dư có phân phối chuẩn.

# 3. Đồ thị ACF của phần dư (ACF of Residuals):
# => Gần như không còn tự tương quan đáng kể (các cột nằm trong dải xanh).

# 4. Đồ thị p-value của kiểm định Ljung-Box:
# => Thấy rằng các giá trị p-value đều > 0.05
# => Nghĩa là phần dư không còn tự tương quan nào.
# => Kết luận: Phần dư giống nhiễu trắng (White Noise).

# --- So sánh AIC của 2 mô hình AR(3) và ARMA(1,1) ---

# Xem các chỉ số của mô hình 2 (giả định là AR(3))
mod2$ICs

# Xem các chỉ số của mô hình 3 (giả định là ARMA(1,1))
mod3$ICs

# Kết luận dựa trên kết quả Console:
# mod3$ICs nhỏ hơn => Chọn mô hình ARMA(1,1)

# --- Dự báo các giá trị tương lai ---
# Cuối cùng, ta có thể dự báo các giá trị tương lai dựa trên mô hình đã chọn
# (Ví dụ dự báo 10 khoảng thời gian tiếp theo)
dubao <- sarima.for(data, n.ahead = 10,p= 1,d = 0,q= 1)

#Các giá trị dự báo 
dubao$pred
#Sai số
dubao$se
# Khoảng tin cậy 95% cho giá trị dự báo
canduoi <- dubao$pred - 1.96 * dubao$se; canduoi
cantren <- dubao$pred + 1.96 * dubao$se; cantren











