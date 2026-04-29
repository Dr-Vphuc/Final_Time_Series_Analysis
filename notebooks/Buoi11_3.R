#11/3/2026
#Qúa trình tuyến tính 
#Qúa trình trung bình trượt bậc q hay còn được gọi là MA(q)


# QUÁ TRÌNH TRUNG BÌNH TRƯỢT

# Quá trình trung bình trượt bậc q: MA(q)
# x_t = w_t + theta_1*w_{t-1} + ... + theta_q*w_{t-q}

# Quá trình trung bình trượt bậc 1: MA(1)
# x_t = w_t + theta*w_{t-1}

# Tính chất:
# Hệ số tương quan khác 0 tại độ trễ 1
# bằng 0 tại các độ trễ lớn hơn 1
# x_t tương quan với x_{t-1} nhưng không tương quan với x_{t-2}, x_{t-3}...

library(TSA)

data(ma1.2.s)
plot(ma1.2.s)
plot(y = ma1.2.s, x = zlag(ma1.2.s,1))
plot(y = ma1.2.s, x = zlag(ma1.2.s, 2))

par(mfrow = c(1, 2))
acf(ma1.2.s)
pacf(ma1.2.s)
# Đối với quá trình MA(1)
# ACF khác 0 tại độ trễ 1, các độ trễ sau bằng 0
# PACF có độ lớn giảm dần

# Quá trình trung bình trượt bậc 2: MA(2)
# x_t = w_t + theta_1*w_{t-1} + theta_2*w_{t-2}

# Tính chất:
# + Hệ số tương quan khác 0 tại độ trễ 1 và 2
#   bằng 0 tại các độ trễ lớn hơn 2
# + x_t tương quan với x_{t-1} và x_{t-2} nhưng không tương quan với x_{t-3}, x_{t-4}...

library(TSA)

data("ma2.s")

par(mfrow = c(2,2))

plot(y = ma2.s, x = zlag(ma2.s,1))

plot(y = ma2.s, x = zlag(ma2.s,2))
plot(y = ma2.s, x = zlag(ma2.s,3))
plot(y = ma2.s, x = zlag(ma2.s,4))

par(mfrow = c(1,2))
acf(ma2.s)
pacf(ma2.s)

# Đối với MA(2)
# ACF khác 0 tại độ trễ 1 và 2, các độ trễ sau bằng 0
# PACF có độ lớn giảm dần

# => Tổng quát: Đối với quá trình MA(q)
# ACF cắt cụt tại độ trễ q 
# (tức là các độ trễ từ 1 đến q thì ACF khác 0
# còn các độ trễ sau q thì ACF bằng 0)

# PACF có độ lớn giảm dần

# Trong quá trình MA thì nếu ta thay \theta bởi 1/(\theta)
# thì giá trị ACF không thay đổi tại các độ trễ

# Điều kiện để quá trình MA khả nghịch là |\theta| < 1
# QUÁ TRÌNH TỰ HỒI QUY

# x_t = phi_1*x_{t-1} + phi_2*x_{t-2} + ... + phi_p*x_{t-p} + w_t

# Quá trình tự hồi quy bậc 1: AR(1)
# x_t = phi*x_{t-1} + w_t

# Tính chất:
# Quá trình AR(1) dừng khi và chỉ khi |phi| < 1

# Đối với các quá trình AR(1) thỏa mãn điều kiện dừng thì
# ACF có độ lớn giảm dần
# PACF khác 0 tại độ trễ 1, còn các độ trễ sau thì bằng 0

data("ar1.2.s")

acf(ar1.2.s); pacf(ar1.2.s)


# Quá trình tự hồi quy bậc 2: AR(2)
# x_t = phi_1*x_{t-1} + phi_2*x_{t-2} + w_t

# Tính chất:
# Quá trình AR(2) dừng khi và chỉ khi nó thỏa mãn các điều kiện sau
# phi_1 + phi_2 < 1
# phi_2 - phi_1 < 1
# |phi_2| < 1

# Đối với các quá trình AR(2) thỏa mãn điều kiện dừng thì
# ACF có độ lớn giảm dần
# PACF khác 0 tại độ trễ 1 và 2, còn các độ trễ sau thì bằng 0

data("ar2.s")

acf(ar2.s)
pacf(ar2.s)


# Tổng quát: Đối với mô hình AR(p) dừng thì
# ACF có độ lớn giảm dần
# PACF cắt cụt tại độ trễ p 
# (tức là p độ trễ đầu tiên khác 0, còn các độ trễ sau bằng 0)

library(tseries)
library(astsa)

data("soi")

plot(soi)

adf.test(soi)

# p-value = 0.01 < 0.05 => Bác bỏ H0, chấp nhận H1
# do đó chuỗi dừng

par(mfrow = c(1,2))

acf(soi)
pacf(soi)

#giáo trình 127/806

