#   25/2/26

#A. Kiểm tra xem chuỗi thời gian có dừng hay không 
#1. adf.test()
library("tseries")
adf.test(oil.price)
#H_0 chuỗi không dừng
#H_1: Chuỗi dừng
#Nếu p_value < 0.05 => Bác bỏ H_0, chấp nhận H_1
                  #tức là chuỗi dừng
#Nếu p_value >= 0.05 =>Chưa đủ cơ sở để bác bỏ H_0, chấp nhận H_1
                  #tức là không chuỗi dừng
#chuối oil.price có p_value = 0.99 > 0.05
#đo đó chuỗi này không dừng

#2.kpss.test()
library("tseries")
kpss.test(oil.price)
#H_0 chuỗi dừng
#H_1: Chuỗi không dừng
#Nếu p_value < 0.05 => Bác bỏ H_0, chấp nhận H_1
#tức là chuỗi không dừng
#Nếu p_value >= 0.05 =>Chưa đủ cơ sở để bác bỏ H_0
#tức là chuỗi dừng
#chuối oil.price có p_value = 0.01 < 0.05 => Bác bỏ H_0
#tức là chuỗi này không dừng


#B. Biến đổi chuỗi thời gian để được chuỗi dừng 

#Cách phổ biến : Sai phân 
#( sử dụng cho các chuỗi cần loại bỏ xu hướng)
diff(oil.price) # sai phân
plot(diff(oil.price))
adf.test(diff(oil.price))
# Thực hiện kiểm định adf với chuỗi oil.price sau khi sai phân 
#thấy rằng p_value = 0.01 < 0.05 => Bác bỏ H0, chấp nhận H_1
#tức là chuỗi sau khi sai phân là chuỗi dừng 

#Cách khác: lấy log
#(sử dụng cho các chuỗi cần ổn định phương sai)


#Hoặc kết hợp giữa 2 phương pháp 
diff(log(oil.price))
adf.test(diff(log(oil.price)))
adf.test(log(diff(oil.price)))#test ngược lại diff -> log

#thực hiện lấy log+diff của chuỗi oil.price
#thấy rằng p_value = 0.01 < 0.05 => Bác bỏ H0, chấp nhận H_1
#tức là chuỗi sau khi lấy log+diff là chuỗi dừng 

#c.Kiểm didngj tính chuẩn của dữ liệu sau khi thực hiện các phương pháp 
#Qqplot
par(mfrow(c(1,3)))
qqnorm(diff(log(oil.price))); qqline(diff(log(oil.price)))
#Kiểm địng Shapiro-Wilk
#H_0 : dữ liệu tuân theo pp chuẩn 
#H_1 : dữ liệu k tuân theo pp chuẩn
#Nếu p_value<0.05 => Bác bỏ H_0, chấp nhận H_1
