#4/3/2026

#Xu hướng được phân loại thành 
#- Xu hướng ngẫu nhiên
#-Xu hướng tất định: đây là xu hướng mà gttb của chuỗi thời gian thay đổi theo 1 quy luật nhất định
#và có thể dự đoán được theo thời gian

#Một số dạng xu hướng tất định phổ biến
#+Xu hướng chu kì/ mùa vụ 
#+Xu hướng tuyến tính: mu_t = beta_0 + beta_1 * t
#+Xu hướng bậc 2, bậc 3 ... mu_t = beta_0 + beta_1*t + beta_2*t^2
#+Xu hướng cosin : mu_t = beta_0 + beta_1* cos(2*pi*f*t) + beta_2*sin(2*pi*f*t)
# f là tần số 
#A. xu hướng tuyến tinh
# A. XU HUONG TUYEN TINH
install.packages("TSA")
library(TSA)
data("rwalk")

rwalk
plot(rwalk)

# 1. Xác định mô hình
# Nhận thấy CTG này có vẻ có xu hướng tuyến tính
# => Ước lượng tham số của mô hình (hồi quy tuyến tính)

n <- length(rwalk)
t <- 1:n
model1 <- lm(rwalk ~ t)
summary(model1)
# => Kết quả

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -1.007888 0.297245  -3.391  0.00126 **
# time(rwalk)  0.134087 0.008475  15.822  < 2e-16 ***

# Ta thấy rằng các tham số của mô hình đều có p-value < 0.05
# tức là các tham số của mô hình đều có ý nghĩa thống kê

# Phương trình:
# rwalk = -1.007888 + 0.134087 * t

# Multiple R-squared: 0.8119, Adjusted R-squared: 0.8086
# F-statistic: 250.3 on 1 and 58 DF, p-value: < 2.2e-16

# Multiple R-squared: 0.8119 cho thấy biến thời gian giải thích được
# 81.19% sự biến thiên của rwalk

# Kiểm định F có p-value < 2.2e-16 < 0.05 cho thấy mô hình hồi quy
# có ý nghĩa thống kê

# 2. Chuẩn đoán mô hình
# Thực hiện phân tích phần dư
#a.kiểm tra bằng đồ thị
# Biểu đồ phần dư theo thời gian
plot(rstudent(model1), type = "l")

# Kiểm tra tính chuẩn bằng QQplot
qqnorm(rstudent(model1))
qqline(rstudent(model1))

# b. Sử dụng các kiểm định thống kê
# Kiểm tra tính chuẩn bằng KD Shapiro-Wilk
shapiro.test(rstudent(model1))

# Ta thấy p-value = 0.8338 > 0.05 => Phần dư của mô hình có PP chuẩn

# Kiểm tra tính độc lập của phần dư bằng kiểm định dãy (Run test)
runs(rstudent(model1))

# Giả thiết H0: Các giá trị trong chuỗi là ngẫu nhiên
# p-value = 0.0266 < 0.05 => bác bỏ H0
# tức là phần dư của mô hình có dấu hiệu không độc lập (phụ thuộc)

# c. Kiểm tra tự tương quan mẫu (ACF)
# để kiểm tra sự phụ thuộc của dữ liệu ở các độ trễ khác nhau
acf(rstudent(model1))

# Ta thấy rằng có một số tự tương quan tại độ trễ VD như 1,2,5,6
# đang vượt quá giới hạn cho phép
# => mô hình chưa xử lí hết được sự phụ thuộc trong dữ liệu

rstudent(model1)

# 3. Dự báo
predict_value <- predict(model1,
                         newdata = data.frame(t = (n + 1):(n + 10)))
predict_value

plot(c(rwalk, predict_value), type = "l")
lines((n + 1):(n + 10), predict_value, col = "red")

# XU HƯỚNG CHU KÌ/MÙA


data("tempdub")
tempdub
plot(tempdub)

# 1. Xác định mô hình
model2 <- lm(tempdub ~ season(tempdub))
summary(model2)

# => Kết quả
# Với mô hình này, R tự động chọn tháng 1 làm gốc và loại nó ra khỏi danh sách.
# Các hệ số đều có p-value < 0.05 tức là đều có ý nghĩa thống kê.

# tempdub = 16.608 + 4.042*I_Feb + 15.867*I_Mar + 29.917*I_Apr + ... + 7.033*I_Dec
# trong đó I ở đây là hàm chỉ số

# Ví dụ: I_Feb = 1 nếu rơi vào tháng 2 và = 0 nếu ngược lại

# Multiple R-squared: 0.9712, Adjusted R-squared: 0.9688
# F-statistic: 405.1 on 11 and 132 DF, p-value: < 2.2e-16

# Multiple R-squared = 0.9712 cho thấy các biến mùa vụ giải thích được
# 97.12% sự biến thiên của tempdub

# Kiểm định F có p-value < 2.2e-16 < 0.05 cho thấy mô hình hồi quy
# có ý nghĩa thống kê

# Ngoài ra chúng ta có thể loại bỏ Intercept trong mô hình hồi quy
model3 <- lm(tempdub ~ season(tempdub) - 1)   # loại bỏ intercept trong mô hình
summary(model3)

# tempdub = 16.608*I_Jan + 20.650*I_Feb + 32.475*I_Mar + 46.525*I_Apr + ... + 23.642*I_Dec
# Trong đó I ở đây là chỉ số
# Như vậy hai mô hình đều cho kết quả dự báo tương tự nhau cho cùng 1 tháng.

# 2. Chuẩn đoán mô hình
# Thực hiện phân tích phần dư

# a. Kiểm tra bằng đồ thị
# Biểu đồ phần dư theo thời gian
plot(rstudent(model2), type = "l")

# Kiểm tra tính chuẩn bằng QQplot
qqnorm(rstudent(model2))
qqline(rstudent(model2))

# b. Sử dụng các kiểm định thống kê
# Kiểm tra tính chuẩn bằng KD Shapiro-Wilk
shapiro.test(rstudent(model2))

# p-value = 0.6954 > 0.05 cho thấy phần dư có phân phối chuẩn

# Kiểm tra tính độc lập của phần dư bằng kiểm định dãy (Run test)
runs(rstudent(model2))

# Giả thiết H0: Các giá trị trong chuỗi là ngẫu nhiên
# p-value = 0.216 > 0.05 => chưa đủ cơ sở để bác bỏ H0
# tức là phần dư của mô hình độc lập

# c. Kiểm tra tự tương quan mẫu (ACF)
# để kiểm tra sự phụ thuộc của dữ liệu ở các độ trễ khác nhau
acf(rstudent(model2))
#Không có tự tương quan tại các độ trễ nào vượt quá giới hạn
#Phần dư của mo hình không có sự phụ thuộc

#phần dư của mô hình giống nhiễu trắng

# XU HƯỚNG COSIN

# 1. Xác định mô hình
model4 <- lm(tempdub ~ harmonic(tempdub, 1))
summary(model4)

# Kết quả
# Ta thấy rằng các hệ số của mô hình đều có p-value < 0.05
# tức là các hệ số đều có ý nghĩa thống kê

# Phương trình:
# tempdub = 46.2660 - 26.7079*cos(2*pi*t) - 2.1697*sin(2*pi*t)

# Multiple R-squared: 0.9639
# Adjusted R-squared: 0.9634

# F-statistic: 1882 on 2 and 141 DF, p-value < 2.2e-16

# Multiple R-squared = 0.9639 cho thấy các biến mùa vụ giải thích được 
# khoảng 96.39% sự biến thiên của tempdub

# Kiểm định F có p-value < 2.2e-16 < 0.05
# => mô hình hồi quy có ý nghĩa thống kê
# 2. Chuẩn đoán mô hình
# Thực hiện phân tích phần dư

# a. Kiểm tra bằng đồ thị
# Biểu đồ phần dư theo thời gian
plot(rstudent(model4), type = "l")

# Kiểm tra tính chuẩn bằng QQplot
qqnorm(rstudent(model4))
qqline(rstudent(model4))

# b. Sử dụng các kiểm định thống kê
# Kiểm tra tính chuẩn bằng KD Shapiro-Wilk
shapiro.test(rstudent(model4))
# p-value = 0.8495 > 0.05 cho thấy phần dư có PP chuẩn

# Kiểm tra tính độc lập của phần dư bằng kiểm định dãy (Run test)
runs(rstudent(model4))
# Giả thiết H0: các giá trị trong chuỗi là ngẫu nhiên
# p-value = 0.112 > 0.05 => chưa đủ cơ sở để bác bỏ H0
# tức là phần dư của mô hình độc lập

# c. Kiểm tra tự tương quan mẫu (ACF)
# để kiểm tra sự phụ thuộc của dữ liệu ở các độ trễ khác nhau
acf(rstudent(model4))
# Có một số tự tương quan tại các độ trễ 2,3,6,9... vượt quá giới hạn
# tức là có dấu hiệu mô hình đang chưa xử lý hết sự phụ thuộc trong DL

# Ta thấy mô hình 2 có phần dư giống nhiễu trắng
# Còn mô hình 4 phần dư vẫn còn sự phụ thuộc
# Do đó nếu xét theo tiêu chí phần dư giống nhiễu trắng thì mô hình 2 tốt hơn

AIC(model2)
AIC(model4)

# Thấy rằng AIC của mô hình 2 nhỏ hơn
# Do đó, ta có thể nói rằng mô hình 2 tốt hơn theo tiêu chí AIC





