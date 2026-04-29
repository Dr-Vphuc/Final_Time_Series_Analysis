#1/4/2026


#1. Phân rã chuỗi thời gian
# Tách CTG thành các thành phần xu hướng, mùa vụ và nhiễu
decompose(co2)
plot(decompose(co2))
install.packages("TSA")
library(TSA)
plot(airpass)
decompose(airpass, type = "multiplicative")
plot(decompose(airpass, type = "multiplicative"))

#2. Hồi quy mùa (đã giới thiệu trong phần hồi quy)
lm(co2 ~ season(co2))

#3. Trung bình trượt mùa vụ
library(forecast)
co2_ma <- ma(co2, order = 12)
co2_ma
plot(co2)
lines(co2_ma, col = "red")

#4. Mô hình SARIMA
plot(co2)

# nhận thấy có tính mùa khá rõ rệt
library(tseries)
adf.test(co2)

# nhận thấy p-value < 0.05 => chuỗi dừng
# => chỉ cần sai phân theo mùa
diff_co2 <- diff(co2, lag = 12)
plot(diff_co2)

#Kiểm tra acf và pacf của chuỗi sau
#sai phân theo mùa với chu kì mùa sa
par(mfrow = c(1,2))
acf(diff(co2, lag = 12), lag.max = 36)
pacf(diff(co2, lag = 12), lag.max = 36)
# Các nhận xét (theo slide)
# ACF giảm dần
# PACF có giá trị đáng kể tại các độ trễ 1, 2, 12 và 24

# => Dự đoán bậc của mô hình
# Với chu kỳ mùa s = 12
# p = 2, d = 0, q = 0
# P = 2, D = 1, Q = 0

# => Mô hình: ARIMA(2,0,0)x(2,1,0)_12

# Fit mô hình SARIMA
library(astsa)  # sarima nằm trong package này
mod1 <- sarima(co2,
               p = 2, d = 0, q = 0,
               P = 2, D = 1, Q = 0,
               S = 12,
               no.constant = TRUE)

# Có thể dự báo bậc của mô hình bằng lệnh auto.arima()
library(forecast)
auto.arima(co2)
mod2 <- sarima(co2,
               p = 1, d = 0, q = 1,
               P = 0, D = 1, Q = 1,S = 12)

mod1$ICs
mod2$ICs
#Dự báo giá trị trong 3 năm tới 
par(mfrow = c(1,1))
sarima.for(co2, n.ahead = 36, p=1,d=0,q=1.)



#10.8
mod2 <- lm(co2 ~ season(co2) + time(co2))
summary(mod2)
acf(residuals(mod2))



