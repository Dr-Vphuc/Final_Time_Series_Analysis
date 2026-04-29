# Vẽ biểu đồ thời gian của dữ liệu.
# Biến đổi dữ liệu.
# + lấy log
# + lấy diff(log)
# + biến đổi lũy thừa
# Xác định bậc của mô hình.
# Ước lượng tham số.
# Chẩn đoán và lựa chọn mô hình

library(astsa)

# mô phỏng 1000 dữ liệu tuần theo mô hình ARIMA(1,1,1)
set.seed(123)
data <- sarima.sim(n = 1000, ar = 0.5, d = 1, ma = 0.2)

plot(data)
# kiểm tra tính dừng
library(tseries)

adf.test(data)
# p-value = 0.3688 > 0.05 => chưa đủ cơ sở bác bỏ H0
# => chuỗi không dừng

adf.test(diff(data))

# p-value = 0.01 < 0.05 => chuỗi dừng
par(mfrow = c(1,2))
acf(diff(data))   # giảm dần
pacf(diff(data))  # gợi ý AR(2)

eacf(diff(data))
# gợi ý mô hình ARMA(1,1) hoặc ARMA(0,3) (đây chính là MA(3))

mod1 <- sarima(diff(data), p = 2, d = 0, q = 0, no.constant = TRUE) #AR(2)
mod2 <- sarima(diff(data), p = 1, d = 0, q = 1, no.constant = TRUE) #AR(2) #ARMA(1,1)
mod3 <- sarima(diff(data), p = 0, d = 0, q = 3, no.constant = TRUE) #MA(3) 
#Phần dư của cả 3 mô hình đều giống nhiễu trắng 
mod1$ICs
mod2$ICs
mod3$ICs










library(TSA)

data("deere3")
plot(deere3)
adf.test(deere3)
#p_value = 0.01 < 0.05 => bác bỏ H_0
#tức chuỗi dừng

acf(deere3)
#acf cắt ở dộ trễ thứ 2
pacf(deere3)
#pacf cắt ở dộ trễ 1

ar1 <- sarima(deere3, p = 1, d = 0, q = 0, no.constant = TRUE) #AR(1)
ar2 <- sarima(deere3, p = 2, d = 0, q = 0, no.constant = TRUE) #AR(2)
#cả 2 mô hình đều khồn có gtri acf của phần dư nào lớn hơn một cách đáng kể các giá trị p_value của kđ Ljung-Box đều > 0.05
#=> trong phần dư không còn tương quan 
#các điểm dữ liệu khá sát đường thẳng chuẩn 
#=> phần dư có tính chuẩn 
#phần dư có tính chất giống nhiễu trắng

ar1$Ics # nhỏ hơn => mô hình AR(1) tốt hơn theo AIC
ar2$ICs
sarima.for(deere3, n.ahead = 10, p =1, d = 0, q = 0)
install.packages("forecast")
library(forecast)
auto.arima(deere3)

