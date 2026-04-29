# BÁO CÁO CHI TIẾT: PHÂN TÍCH CHUỖI THỜI GIAN

> Tổng hợp lý thuyết và minh họa code R từ 6 buổi học trong [notebooks/](notebooks/).
> Tác giả: tổng hợp từ bài giảng — 2026.

---

## Mục lục

- [Phần I — Khái niệm cơ bản & Tiền xử lý chuỗi](#phần-i--khái-niệm-cơ-bản--tiền-xử-lý-chuỗi)
  - [Chương 1. Tính dừng của chuỗi thời gian](#chương-1-tính-dừng-của-chuỗi-thời-gian)
  - [Chương 2. Biến đổi chuỗi để đạt tính dừng](#chương-2-biến-đổi-chuỗi-để-đạt-tính-dừng)
- [Phần II — Mô hình xu hướng tất định](#phần-ii--mô-hình-xu-hướng-tất-định)
  - [Chương 3. Phân loại xu hướng](#chương-3-phân-loại-xu-hướng)
  - [Chương 4. Hồi quy xu hướng tuyến tính](#chương-4-hồi-quy-xu-hướng-tuyến-tính)
  - [Chương 5. Hồi quy xu hướng mùa vụ](#chương-5-hồi-quy-xu-hướng-mùa-vụ)
  - [Chương 6. Hồi quy xu hướng cosin (điều hòa)](#chương-6-hồi-quy-xu-hướng-cosin-điều-hòa)
- [Phần III — Quá trình ngẫu nhiên tuyến tính](#phần-iii--quá-trình-ngẫu-nhiên-tuyến-tính)
  - [Chương 7. Quá trình trung bình trượt MA(q)](#chương-7-quá-trình-trung-bình-trượt-maq)
  - [Chương 8. Quá trình tự hồi quy AR(p)](#chương-8-quá-trình-tự-hồi-quy-arp)
  - [Chương 9. Mô hình ARMA(p, q)](#chương-9-mô-hình-armap-q)
- [Phần IV — Mô hình ARIMA & SARIMA](#phần-iv--mô-hình-arima--sarima)
  - [Chương 10. Mô hình ARIMA(p, d, q)](#chương-10-mô-hình-arimap-d-q)
  - [Chương 11. Phân rã chuỗi thời gian](#chương-11-phân-rã-chuỗi-thời-gian)
  - [Chương 12. Mô hình SARIMA(p, d, q)(P, D, Q)_s](#chương-12-mô-hình-sarimap-d-qp-d-q_s)
- [Phần V — Phụ lục](#phần-v--phụ-lục)
  - [A. Tổng hợp các package R sử dụng](#a-tổng-hợp-các-package-r-sử-dụng)
  - [B. Tổng hợp các kiểm định thống kê](#b-tổng-hợp-các-kiểm-định-thống-kê)
  - [C. Bảng nhận diện ACF/PACF cho MA / AR / ARMA](#c-bảng-nhận-diện-acfpacf-cho-ma--ar--arma)
  - [D. Sơ đồ quy trình Box-Jenkins](#d-sơ-đồ-quy-trình-box-jenkins)
  - [E. Danh mục datasets sử dụng](#e-danh-mục-datasets-sử-dụng)

---

# Phần I — Khái niệm cơ bản & Tiền xử lý chuỗi

> Nguồn: [Buoi25_2.R](notebooks/Buoi25_2.R)

## Chương 1. Tính dừng của chuỗi thời gian

### 1.1. Định nghĩa chuỗi dừng

Một chuỗi thời gian $\{X_t\}$ được gọi là **dừng yếu (weakly stationary)** nếu:

1. $E[X_t] = \mu$ — kỳ vọng không đổi theo thời gian.
2. $\operatorname{Var}(X_t) = \sigma^2 < \infty$ — phương sai hữu hạn, không đổi.
3. $\operatorname{Cov}(X_t, X_{t+h}) = \gamma(h)$ — hiệp phương sai chỉ phụ thuộc khoảng trễ $h$, không phụ thuộc $t$.

Chuỗi **dừng mạnh (strictly stationary)** khi phân phối đồng thời của $(X_{t_1}, \ldots, X_{t_k})$ bằng với $(X_{t_1+h}, \ldots, X_{t_k+h})$ với mọi $h$.

Tính dừng là **điều kiện tiên quyết** để áp dụng các mô hình AR, MA, ARMA, ARIMA.

### 1.2. Kiểm định Augmented Dickey-Fuller (ADF)

Mô hình ADF kiểm tra sự tồn tại của **nghiệm đơn vị (unit root)** trong phương trình hồi quy:

$$\Delta X_t = \alpha + \beta t + \gamma X_{t-1} + \sum_{i=1}^{p} \delta_i \, \Delta X_{t-i} + \varepsilon_t$$

- $H_0$: $\gamma = 0$ — chuỗi **không dừng** (có nghiệm đơn vị).
- $H_1$: $\gamma < 0$ — chuỗi **dừng**.

**Quy tắc quyết định**:
- Nếu $p\text{-value} < 0.05$ → bác bỏ $H_0$ → chuỗi **dừng**.
- Nếu $p\text{-value} \geq 0.05$ → chưa đủ cơ sở bác bỏ $H_0$ → chuỗi **không dừng**.

```r
library("tseries")
adf.test(oil.price)
```

**Kết quả minh họa**: với chuỗi `oil.price`, $p\text{-value} = 0.99 > 0.05$ → chuỗi không dừng.

### 1.3. Kiểm định KPSS (Kwiatkowski–Phillips–Schmidt–Shin)

KPSS có giả thuyết **ngược lại** ADF:

- $H_0$: chuỗi **dừng** (xung quanh một mức hoặc xu hướng tất định).
- $H_1$: chuỗi **không dừng**.

**Quy tắc quyết định**:
- Nếu $p\text{-value} < 0.05$ → bác bỏ $H_0$ → chuỗi **không dừng**.
- Nếu $p\text{-value} \geq 0.05$ → chưa đủ cơ sở bác bỏ $H_0$ → chuỗi **dừng**.

```r
library("tseries")
kpss.test(oil.price)
```

**Kết quả minh họa**: với chuỗi `oil.price`, $p\text{-value} = 0.01 < 0.05$ → bác bỏ $H_0$ → chuỗi không dừng (cùng kết luận với ADF).

### 1.4. So sánh ADF vs KPSS

| Tiêu chí | ADF | KPSS |
|----------|-----|------|
| $H_0$ | Không dừng | Dừng |
| $H_1$ | Dừng | Không dừng |
| Bác bỏ $H_0$ ⇒ | **Dừng** | **Không dừng** |

Trong thực hành, **dùng cả hai** và đối chiếu để có kết luận đáng tin cậy nhất.

> 📌 **Ý nghĩa & Khuyến nghị sử dụng (bổ sung)**
>
> - **Tại sao cần kiểm tra tính dừng?** Hầu hết các mô hình kinh điển (AR, MA, ARMA, ARIMA) đều giả định chuỗi dừng. Nếu áp dụng trên chuỗi không dừng, ước lượng tham số sẽ chệch và dự báo trở nên không tin cậy ("hồi quy giả tạo" — spurious regression).
> - **ADF** kiểm tra sự tồn tại của *unit root* ($\gamma = 0$). Nó **mạnh** khi mẫu lớn, nhưng **yếu** (low power) khi mẫu nhỏ hoặc khi chuỗi gần giáp ranh dừng/không dừng.
> - **KPSS** dựa trên mô hình phân rã chuỗi thành xu hướng tất định + bước ngẫu nhiên + sai số dừng. Ưu thế lớn của KPSS là **phát hiện được trend-stationarity** (chuỗi dừng quanh xu hướng tất định), điều mà ADF có thể bỏ sót.
> - **Quy tắc kết hợp 4 trường hợp**:
>   - ADF: dừng & KPSS: dừng → Chuỗi **đã dừng**, dùng được ngay.
>   - ADF: không dừng & KPSS: không dừng → Chuỗi **không dừng**, cần sai phân.
>   - KPSS bác bỏ, ADF không bác bỏ → **Difference-stationary** → cần `diff()`.
>   - KPSS không bác bỏ, ADF bác bỏ → **Trend-stationary** → cần loại bỏ xu hướng tất định bằng hồi quy.
> - **Khi nào KHÔNG nên** chỉ dùng một kiểm định? Khi mẫu nhỏ ($n < 100$) hoặc chuỗi có thay đổi cấu trúc (structural break) — nên kết hợp với kiểm định Phillips-Perron, Zivot-Andrews.

---

## Chương 2. Biến đổi chuỗi để đạt tính dừng

### 2.1. Sai phân (differencing)

**Sai phân bậc 1** loại bỏ xu hướng tuyến tính:

$$\nabla X_t = X_t - X_{t-1}$$

**Sai phân bậc $d$**: áp dụng phép sai phân $d$ lần liên tiếp.

```r
diff(oil.price)            # sai phân bậc 1
plot(diff(oil.price))
adf.test(diff(oil.price))
# p-value = 0.01 < 0.05 => chuỗi sau sai phân là dừng
```

### 2.2. Lấy logarit

Phép biến đổi $\log(X_t)$ giúp **ổn định phương sai** (đặc biệt với các chuỗi có phương sai tăng theo mức giá trị, ví dụ chuỗi tài chính, lượng hành khách).

### 2.3. Kết hợp `diff(log(.))`

Đây là biến đổi rất phổ biến: vừa ổn định phương sai (qua `log`), vừa loại bỏ xu hướng (qua `diff`). Đại lượng `diff(log(X_t))` xấp xỉ **tỷ lệ tăng trưởng** (return) của chuỗi.

```r
diff(log(oil.price))
adf.test(diff(log(oil.price)))
# p-value = 0.01 < 0.05 => chuỗi sau biến đổi log+diff là dừng
```

> **Lưu ý**: thứ tự `diff(log())` và `log(diff())` cho kết quả khác nhau. Thứ tự đúng/được khuyến nghị là `diff(log(.))` vì `diff(.)` có thể tạo giá trị âm khiến `log(.)` không xác định.

### 2.4. Kiểm tra tính chuẩn của chuỗi sau biến đổi

#### a) QQ-plot
```r
qqnorm(diff(log(oil.price)))
qqline(diff(log(oil.price)))
```
Nếu các điểm nằm sát đường thẳng tham chiếu → dữ liệu xấp xỉ phân phối chuẩn.

#### b) Kiểm định Shapiro-Wilk

- $H_0$: dữ liệu tuân theo phân phối chuẩn.
- $H_1$: dữ liệu không tuân theo phân phối chuẩn.

Nếu $p\text{-value} < 0.05$ → bác bỏ $H_0$ → dữ liệu **không** tuân phân phối chuẩn.

```r
shapiro.test(diff(log(oil.price)))
```

> 📌 **Ý nghĩa & Khuyến nghị sử dụng (bổ sung)**
>
> **Sai phân `diff()`**
> - **Khái niệm**: phép biến đổi $\nabla X_t = X_t - X_{t-1}$ giúp loại bỏ xu hướng và mùa vụ ở mức tuyến tính, đưa chuỗi không dừng (do *xu hướng ngẫu nhiên*) về chuỗi dừng.
> - **Ưu điểm**: đơn giản, không cần ước lượng tham số xu hướng; phù hợp với dữ liệu kinh tế-tài chính (giá cổ phiếu, tỷ giá), chuỗi có **xu hướng ngẫu nhiên** (random walk).
> - **Nhược điểm**: mất một quan sát đầu tiên mỗi lần sai phân; không tách biệt được thành phần xu hướng và thành phần ngẫu nhiên; có thể sai phân *quá mức* (over-differencing) → tạo ra tự tương quan âm giả.
> - **Khi nào nên dùng**: chuỗi có *unit root* (ADF không bác bỏ $H_0$), chuỗi có xu hướng nhưng không phải tất định.
> - **Khi nào KHÔNG nên dùng**: chuỗi đã dừng — sai phân không cần thiết sẽ làm tăng phương sai phần dư; chuỗi chỉ có *trend-stationary* — nên dùng hồi quy loại trend thay vì sai phân.
>
> **Lấy logarit `log()`**
> - **Khái niệm**: biến đổi đơn điệu giúp **ổn định phương sai** khi phương sai tăng theo mức của chuỗi (heteroscedasticity nhân tính).
> - **Ý nghĩa**: $\nabla \log(X_t) \approx \dfrac{X_t - X_{t-1}}{X_{t-1}}$ — chính là **tỷ lệ tăng trưởng** (growth rate / log-return). Hệ số góc của đường xu hướng trên log-scale = tốc độ tăng trưởng trung bình hàng kỳ.
> - **Ưu điểm**: giảm độ lệch (skewness), tuyến tính hóa tăng trưởng theo cấp số nhân, dễ diễn giải.
> - **Nhược điểm**: **không áp dụng được với giá trị âm hoặc 0**; nếu phương sai vốn đã ổn định, lấy log có thể làm xấu mô hình.
> - **Khi nào nên dùng**: chuỗi tài chính, dân số, sản lượng, giá hàng hóa — bất cứ khi nào dao động tỉ lệ với mức.
> - **Khi nào KHÔNG nên dùng**: chuỗi có giá trị 0/âm; chuỗi có phương sai ổn định; chuỗi mà người dùng cần dự báo trên thang đo gốc một cách chính xác (vì biến đổi ngược $\exp$ tạo bias).
>
> **Kết hợp `diff(log(.))`**: cách phổ biến nhất với chuỗi tài chính — vừa ổn định phương sai vừa loại xu hướng. Kết quả gần đúng = chuỗi log-return.

---

# Phần II — Mô hình xu hướng tất định

> Nguồn: [Buoi4_3.R](notebooks/Buoi4_3.R)

## Chương 3. Phân loại xu hướng

Xu hướng (trend) trong chuỗi thời gian được chia làm hai loại:

- **Xu hướng ngẫu nhiên (stochastic trend)**: do tác động của nhiễu tích lũy theo thời gian (ví dụ: random walk).
- **Xu hướng tất định (deterministic trend)**: giá trị trung bình của chuỗi thay đổi theo một quy luật xác định và **dự đoán được** theo thời gian.

### Các dạng xu hướng tất định phổ biến

| Dạng | Phương trình |
|------|--------------|
| Tuyến tính | $\mu_t = \beta_0 + \beta_1 t$ |
| Bậc 2 / bậc 3 | $\mu_t = \beta_0 + \beta_1 t + \beta_2 t^2 + \cdots$ |
| Mùa vụ (chỉ số) | $\mu_t = \sum_{j=1}^{s} \beta_j \mathbb{I}_{\{t \in \text{tháng}_j\}}$ |
| Cosin (điều hòa) | $\mu_t = \beta_0 + \beta_1 \cos(2\pi f t) + \beta_2 \sin(2\pi f t)$ |

> 📌 **Ý nghĩa & Phân biệt sâu (bổ sung)**
>
> - **Xu hướng tất định** ngầm giả định *độ dốc xu hướng không đổi theo thời gian* — tốc độ tăng trưởng trong tương lai vẫn bằng tốc độ tăng trưởng quá khứ. Ưu điểm: mô hình đơn giản, khoảng tin cậy hẹp, dễ diễn giải.
> - **Xu hướng ngẫu nhiên** cho phép độ dốc *thay đổi theo thời gian* — tăng trưởng ước lượng chỉ là trung bình lịch sử, không nhất thiết duy trì. Nhược điểm: khoảng tin cậy dự báo *rộng hơn nhiều* vì sai số không dừng; nhưng phù hợp hơn với hiện thực của hầu hết chuỗi kinh tế.
> - **Vấn đề thực tế**: với mẫu hữu hạn, gần như luôn có một mô hình xu hướng tất định và một mô hình xu hướng ngẫu nhiên *fit data tương đương nhau* → cần kết hợp ADF/KPSS, kiểm tra phần dư, và **kiến thức chuyên ngành** để chọn loại trend.
> - **Khi nào nên dùng xu hướng tất định**: chuỗi vật lý/kỹ thuật có cơ sở lý thuyết về xu hướng (nhiệt độ theo mùa, chu kỳ thiên văn), chuỗi ngắn hạn ổn định.
> - **Khi nào KHÔNG nên dùng**: chuỗi tài chính, kinh tế dài hạn — nên dùng ARIMA (xu hướng ngẫu nhiên qua sai phân) thay vì hồi quy theo $t$.

---

## Chương 4. Hồi quy xu hướng tuyến tính

### 4.1. Mô hình

$$X_t = \beta_0 + \beta_1 t + e_t, \qquad t = 1, 2, \ldots, n$$

trong đó $e_t$ là sai số ngẫu nhiên (kỳ vọng giả định là nhiễu trắng).

### 4.2. Ước lượng tham số

Sử dụng dữ liệu `rwalk` từ package `TSA`:

```r
library(TSA)
data("rwalk")
plot(rwalk)

n <- length(rwalk)
t <- 1:n
model1 <- lm(rwalk ~ t)
summary(model1)
```

### 4.3. Đọc kết quả `summary()`

```
Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) -1.007888 0.297245  -3.391  0.00126 **
t            0.134087 0.008475  15.822  < 2e-16 ***

Multiple R-squared: 0.8119,  Adjusted R-squared: 0.8086
F-statistic: 250.3 on 1 and 58 DF,  p-value: < 2.2e-16
```

- Cả hai hệ số đều có $p\text{-value} < 0.05$ → có ý nghĩa thống kê.
- Phương trình ước lượng: $\widehat{X_t} = -1.0079 + 0.1341 \, t$.
- $R^2 = 0.8119$ → biến thời gian giải thích **81.19%** sự biến thiên của `rwalk`.
- Kiểm định $F$ có $p\text{-value} < 2.2 \times 10^{-16}$ → mô hình có ý nghĩa thống kê.

### 4.4. Chẩn đoán phần dư

#### a) Kiểm tra bằng đồ thị

```r
plot(rstudent(model1), type = "l")   # phần dư chuẩn hóa theo thời gian
qqnorm(rstudent(model1))
qqline(rstudent(model1))
```

#### b) Kiểm định tính chuẩn — Shapiro-Wilk

```r
shapiro.test(rstudent(model1))
# p-value = 0.8338 > 0.05 => phần dư có phân phối chuẩn
```

#### c) Kiểm định tính độc lập — Run test

- $H_0$: các giá trị trong chuỗi là **ngẫu nhiên** (độc lập).
- $H_1$: các giá trị có sự phụ thuộc.

```r
runs(rstudent(model1))
# p-value = 0.0266 < 0.05 => bác bỏ H0 => phần dư có dấu hiệu phụ thuộc
```

#### d) ACF của phần dư

```r
acf(rstudent(model1))
```
Nếu nhiều giá trị ACF vượt quá dải tin cậy (vạch xanh) → mô hình **chưa xử lý hết** sự phụ thuộc trong dữ liệu (cần mô hình phức tạp hơn).

### 4.5. Dự báo

```r
predict_value <- predict(model1, newdata = data.frame(t = (n + 1):(n + 10)))
plot(c(rwalk, predict_value), type = "l")
lines((n + 1):(n + 10), predict_value, col = "red")
```

> 📌 **Hồi quy tuyến tính theo thời gian — Ưu / Nhược (bổ sung)**
>
> - **Ý nghĩa**: mô hình đơn giản nhất biểu diễn xu hướng — phù hợp cho nghiên cứu mô tả ban đầu.
> - **Ưu điểm**: dễ ước lượng (OLS), dễ diễn giải, nhanh; cho khoảng tin cậy hẹp.
> - **Nhược điểm**:
>   - Giả định độ dốc cố định trên toàn chuỗi → sai khi xu hướng thay đổi.
>   - Phần dư của mô hình thường còn tự tương quan (ví dụ với `rwalk` trong code, run-test cho $p = 0.027 < 0.05$) → cần thêm tầng AR/MA cho phần dư.
>   - Sai khi áp dụng cho chuỗi có *unit root* (kết quả $R^2$ cao một cách giả tạo — spurious regression).
> - **Khi nào nên dùng**: chuỗi rõ ràng tăng/giảm tuyến tính, không có mùa vụ, mẫu vừa phải, cần dự báo ngắn hạn.
> - **Khi nào KHÔNG nên dùng**: chuỗi có *unit root* (luôn kiểm tra ADF/KPSS trước!), chuỗi có mùa vụ rõ rệt, dự báo dài hạn.

---

## Chương 5. Hồi quy xu hướng mùa vụ

### 5.1. Mô hình biến chỉ số

Với dữ liệu nhiệt độ tháng `tempdub` (chu kỳ $s = 12$):

$$X_t = \beta_0 + \sum_{j=2}^{12} \beta_j \, \mathbb{I}_{\{\text{tháng}_j\}} + e_t$$

trong đó $\mathbb{I}_{\{\text{tháng}_j\}} = 1$ nếu $t$ rơi vào tháng $j$, ngược lại bằng 0. R **mặc định** chọn tháng 1 làm gốc và bỏ khỏi danh sách (để tránh đa cộng tuyến).

```r
data("tempdub")
plot(tempdub)
model2 <- lm(tempdub ~ season(tempdub))
summary(model2)
```

Kết quả: $R^2 = 0.9712$, $F$-statistic $= 405.1$, $p\text{-value} < 2.2 \times 10^{-16}$.
Mô hình giải thích **97.12%** biến thiên của nhiệt độ.

### 5.2. Loại bỏ intercept

```r
model3 <- lm(tempdub ~ season(tempdub) - 1)
summary(model3)
```
Cú pháp `- 1` loại bỏ hằng số chặn → mỗi tháng có một hệ số riêng (nhiệt độ trung bình tháng đó).

Hai mô hình tương đương về dự báo: cùng cho ra trung bình mẫu của tháng tương ứng.

### 5.3. Chẩn đoán mô hình

```r
plot(rstudent(model2), type = "l")
qqnorm(rstudent(model2)); qqline(rstudent(model2))
shapiro.test(rstudent(model2))   # p-value = 0.6954 > 0.05 => phần dư chuẩn
runs(rstudent(model2))           # p-value = 0.216 > 0.05  => phần dư độc lập
acf(rstudent(model2))            # không có ACF vượt giới hạn
```

→ Phần dư của mô hình mùa vụ **giống nhiễu trắng**: có chuẩn, có độc lập, không tự tương quan.

> 📌 **Hồi quy mùa vụ chỉ số — Ưu / Nhược (bổ sung)**
>
> - **Ý nghĩa**: ước lượng riêng một hệ số cho mỗi mùa (ví dụ 12 hệ số cho dữ liệu tháng) → mô hình rất linh hoạt, "thuộc lòng" pattern mùa.
> - **Ưu điểm**: cực kỳ linh hoạt — bắt được mọi hình dạng mùa vụ kể cả các đỉnh nhọn / không đối xứng; dễ diễn giải (mỗi hệ số = giá trị trung bình của mùa đó).
> - **Nhược điểm**:
>   - Cần ước lượng nhiều tham số (s tham số) → tốn bậc tự do, kém ổn định khi chuỗi ngắn.
>   - Không hiệu quả với mùa vụ chu kỳ dài (ví dụ $s = 365$ cho dữ liệu ngày-trong-năm) → khi đó nên dùng hồi quy điều hòa (Fourier).
>   - Giả định mùa vụ *cố định* qua các năm — không bắt được mùa vụ thay đổi dần.
> - **Khi nào nên dùng**: chuỗi có chu kỳ ngắn ($s \leq 24$), mùa vụ ổn định, mẫu nhiều năm.
> - **Khi nào KHÔNG nên dùng**: chuỗi rất ngắn (< 2 chu kỳ); chu kỳ rất dài; mùa vụ đang biến đổi (cần mô hình động như TBATS, Prophet).

---

## Chương 6. Hồi quy xu hướng cosin (điều hòa)

### 6.1. Mô hình điều hòa

$$X_t = \beta_0 + \beta_1 \cos(2\pi f t) + \beta_2 \sin(2\pi f t) + e_t$$

trong đó $f$ là tần số. Với dữ liệu tháng có chu kỳ năm, $f = 1$ (đơn vị thời gian được chuẩn hóa theo năm).

### 6.2. Ước lượng

```r
model4 <- lm(tempdub ~ harmonic(tempdub, 1))
summary(model4)
```

Phương trình ước lượng:
$$\widehat{X_t} = 46.266 - 26.708 \cos(2\pi t) - 2.170 \sin(2\pi t)$$

$R^2 = 0.9639$, $F\text{-statistic} = 1882$, $p\text{-value} < 2.2 \times 10^{-16}$.

### 6.3. So sánh các mô hình bằng AIC

```r
AIC(model2)   # mô hình mùa vụ chỉ số
AIC(model4)   # mô hình cosin
```

**Tiêu chí AIC** (Akaike Information Criterion):

$$\text{AIC} = -2 \log L + 2k$$

trong đó $L$ là hàm hợp lý cực đại, $k$ là số tham số. **AIC nhỏ hơn → mô hình tốt hơn**.

Trong ví dụ này:
- AIC của `model2` (mùa vụ chỉ số) **nhỏ hơn** AIC của `model4` (cosin).
- ACF phần dư của `model2` không vượt giới hạn, còn `model4` vẫn còn tự tương quan tại các độ trễ 2, 3, 6, 9...

→ **Kết luận**: `model2` tốt hơn cả về AIC lẫn tính nhiễu trắng của phần dư.

> 📌 **Hồi quy điều hòa (Fourier) — Ưu / Nhược (bổ sung)**
>
> - **Ý nghĩa**: dùng tổ hợp $\sin / \cos$ để xấp xỉ hàm tuần hoàn — bậc càng cao càng linh hoạt, nhưng cũng càng dễ overfit.
> - **Ưu điểm**:
>   - **Tiết kiệm tham số** với chu kỳ dài: chỉ cần $2K$ tham số ($K$ cặp Fourier) thay vì $s - 1$ tham số mùa vụ chỉ số.
>   - **Hiệu quả tính toán** với chu kỳ rất lớn (ngày trong năm $s = 365$, nửa giờ trong tuần $s = 336$).
>   - Cho phép **đồng thời nhiều chu kỳ** (ví dụ điện năng có chu kỳ ngày + tuần + năm).
> - **Nhược điểm**:
>   - Giả định mùa vụ *đối xứng và mượt* — kém với pattern có **đỉnh nhọn**.
>   - Bậc Fourier $K$ là siêu tham số khó chọn (quá nhỏ → underfit, quá lớn → overfit).
>   - Mùa vụ giả định *cố định theo thời gian*.
> - **Khi nào nên dùng**: chu kỳ dài, nhiều chu kỳ chồng chất, dữ liệu tần suất cao.
> - **Khi nào KHÔNG nên dùng**: mùa vụ có dạng nhọn / bất đối xứng (ví dụ doanh số dịp lễ); chuỗi ngắn; mùa vụ thay đổi nhiều theo thời gian.

> 📌 **Tiêu chí AIC (bổ sung)**
>
> - **AIC** = $-2 \log L + 2k$. **BIC** = $-2 \log L + k \log n$.
> - AIC ưu tiên *khả năng dự báo*; BIC ưu tiên *tính tiết kiệm*. BIC phạt nặng hơn → thường chọn mô hình ít tham số hơn AIC.
> - **Cảnh báo**: AIC chỉ so sánh được giữa các mô hình *fit trên cùng dữ liệu* và cùng phương pháp ước lượng. Không so sánh được giữa mô hình ARIMA(p,d,q) với khác $d$ (vì khác cỡ mẫu sau sai phân).

---

# Phần III — Quá trình ngẫu nhiên tuyến tính

## Chương 7. Quá trình trung bình trượt MA(q)

> Nguồn: [Buoi11_3.R](notebooks/Buoi11_3.R)

### 7.1. Định nghĩa

Quá trình **trung bình trượt bậc $q$** — viết tắt $\mathrm{MA}(q)$:

$$X_t = w_t + \theta_1 w_{t-1} + \theta_2 w_{t-2} + \cdots + \theta_q w_{t-q}$$

trong đó $\{w_t\} \overset{\text{iid}}{\sim} \mathcal{N}(0, \sigma_w^2)$ là **nhiễu trắng** (white noise).

### 7.2. MA(1)

$$X_t = w_t + \theta w_{t-1}$$

**Tính chất ACF**:

$$\rho(1) = \frac{\theta}{1 + \theta^2}, \qquad \rho(k) = 0 \;\; \forall k \geq 2$$

→ ACF chỉ khác 0 tại độ trễ 1, các độ trễ sau bằng 0.
→ PACF có độ lớn **giảm dần**.

```r
library(TSA)
data(ma1.2.s)
plot(ma1.2.s)
plot(y = ma1.2.s, x = zlag(ma1.2.s, 1))   # tương quan với x_{t-1}
plot(y = ma1.2.s, x = zlag(ma1.2.s, 2))   # gần như không tương quan với x_{t-2}

par(mfrow = c(1, 2))
acf(ma1.2.s)
pacf(ma1.2.s)
```

### 7.3. MA(2)

$$X_t = w_t + \theta_1 w_{t-1} + \theta_2 w_{t-2}$$

→ ACF khác 0 tại độ trễ 1 và 2, **bằng 0 từ độ trễ 3 trở đi**.
→ PACF có độ lớn giảm dần.

```r
data("ma2.s")
par(mfrow = c(2, 2))
plot(y = ma2.s, x = zlag(ma2.s, 1))
plot(y = ma2.s, x = zlag(ma2.s, 2))
plot(y = ma2.s, x = zlag(ma2.s, 3))
plot(y = ma2.s, x = zlag(ma2.s, 4))

par(mfrow = c(1, 2))
acf(ma2.s)
pacf(ma2.s)
```

### 7.4. Quy tắc tổng quát MA(q)

| Hàm | Hành vi |
|-----|---------|
| ACF | **Cắt cụt tại độ trễ $q$** (khác 0 từ 1 đến $q$, bằng 0 từ $q+1$) |
| PACF | Giảm dần (tail-off) |

### 7.5. Tính khả nghịch (invertibility)

Với MA(1), nếu thay $\theta$ bằng $1/\theta$ thì ACF không thay đổi → mô hình **không nhận dạng được duy nhất**. Để có nghiệm duy nhất, ta yêu cầu điều kiện khả nghịch:

$$|\theta| < 1$$

Tổng quát hơn, $\mathrm{MA}(q)$ khả nghịch khi và chỉ khi **mọi nghiệm** của đa thức MA $\theta(z) = 1 + \theta_1 z + \cdots + \theta_q z^q$ đều nằm **ngoài** đường tròn đơn vị.

> 📌 **MA(q) — Ý nghĩa & Khi nào nên dùng (bổ sung)**
>
> - **Ý nghĩa**: $X_t$ là tổ hợp tuyến tính của các *cú sốc* (innovations / shocks) gần đây. Cú sốc tại $t-q-1$ trở về trước **không** ảnh hưởng đến $X_t$ → MA có "trí nhớ ngắn".
> - **Ưu điểm**:
>   - **Luôn dừng** với mọi giá trị tham số — không cần kiểm tra điều kiện dừng.
>   - Dễ chẩn đoán bậc qua ACF (cắt cụt rõ rệt tại $q$).
>   - Phù hợp với chuỗi mà ảnh hưởng của thông tin/sự kiện *tan biến nhanh* sau vài bước (ví dụ sai số đo lường, hiệu ứng tin tức ngắn hạn).
> - **Nhược điểm**:
>   - Cần điều kiện *khả nghịch* để mô hình duy nhất.
>   - Sau $q$ bước dự báo, MA thuần trở về dự báo bằng mức trung bình → **không phù hợp dự báo dài hạn**.
>   - Không bắt được phụ thuộc dài hạn.
> - **Khi nào nên dùng**: chuỗi dừng, ACF cắt cụt rõ ràng, dự báo ngắn hạn (≤ q bước).
> - **Khi nào KHÔNG nên dùng**: dự báo nhiều bước về tương lai; chuỗi có phụ thuộc dài hạn; chuỗi không dừng (cần ARIMA).

---

## Chương 8. Quá trình tự hồi quy AR(p)

### 8.1. Định nghĩa

Quá trình **tự hồi quy bậc $p$** — viết tắt $\mathrm{AR}(p)$:

$$X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + \cdots + \phi_p X_{t-p} + w_t$$

### 8.2. AR(1)

$$X_t = \phi X_{t-1} + w_t$$

**Điều kiện dừng**: $|\phi| < 1$.

**Tính chất**:
- ACF: $\rho(k) = \phi^k$ → có độ lớn **giảm dần** (theo cấp số nhân).
- PACF: khác 0 tại độ trễ 1, bằng 0 với mọi độ trễ $\geq 2$.

```r
data("ar1.2.s")
acf(ar1.2.s); pacf(ar1.2.s)
```

### 8.3. AR(2)

$$X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + w_t$$

**Điều kiện dừng** (3 bất đẳng thức):

$$\phi_1 + \phi_2 < 1, \qquad \phi_2 - \phi_1 < 1, \qquad |\phi_2| < 1$$

→ ACF giảm dần.
→ PACF khác 0 tại độ trễ 1 và 2, bằng 0 từ độ trễ 3 trở đi.

```r
data("ar2.s")
acf(ar2.s)
pacf(ar2.s)
```

### 8.4. Quy tắc tổng quát AR(p)

| Hàm | Hành vi |
|-----|---------|
| ACF | Giảm dần (tail-off) |
| PACF | **Cắt cụt tại độ trễ $p$** (khác 0 từ 1 đến $p$, bằng 0 từ $p+1$) |

Tổng quát, $\mathrm{AR}(p)$ dừng khi và chỉ khi mọi nghiệm của đa thức AR $\phi(z) = 1 - \phi_1 z - \cdots - \phi_p z^p$ đều nằm **ngoài** đường tròn đơn vị.

### 8.5. Ví dụ thực tế: chuỗi `soi`

```r
library(tseries)
library(astsa)
data("soi")
plot(soi)

adf.test(soi)
# p-value = 0.01 < 0.05 => chuỗi dừng

par(mfrow = c(1, 2))
acf(soi)
pacf(soi)
```

→ Quan sát ACF/PACF để chọn bậc của mô hình AR phù hợp.

> 📌 **AR(p) — Ý nghĩa & Khi nào nên dùng (bổ sung)**
>
> - **Ý nghĩa**: $X_t$ là tổ hợp tuyến tính của *các giá trị quá khứ của chính nó*. Phụ thuộc dài hạn được mã hóa qua các hệ số $\phi_i$.
> - **Ưu điểm**:
>   - Diễn giải trực quan: "$X_t$ phụ thuộc vào $X_{t-1}, \ldots, X_{t-p}$".
>   - Dễ dự báo nhiều bước (đệ quy theo công thức).
>   - Bắt được phụ thuộc dài hạn dạng giảm theo cấp số nhân.
>   - Có thể mô tả chuỗi có **chu kỳ giả** (pseudo-cycles) khi nghiệm đa thức $\phi(z)$ là số phức.
> - **Nhược điểm**:
>   - Phải thỏa mãn điều kiện dừng (ràng buộc tham số).
>   - Không biểu diễn được phần "cú sốc" có ảnh hưởng cấu trúc lên $X_t$.
>   - Số tham số tăng tuyến tính theo $p$ → dễ overfit nếu chọn $p$ quá lớn.
> - **Khi nào nên dùng**: chuỗi dừng có **PACF cắt cụt** rõ ràng; chuỗi kinh tế/khí tượng có quán tính.
> - **Khi nào KHÔNG nên dùng**: chuỗi không dừng; PACF không cắt cụt rõ ràng (nên thử ARMA); chuỗi có thay đổi chế độ (regime change) — AR tuyến tính sẽ thất bại.

---

## Chương 9. Mô hình ARMA(p, q)

> Nguồn: [Buoi18_3.R](notebooks/Buoi18_3.R)

### 9.1. Định nghĩa

Quá trình **trung bình trượt tự hồi quy** $\mathrm{ARMA}(p, q)$ kết hợp cả AR và MA:

$$X_t = \phi_1 X_{t-1} + \cdots + \phi_p X_{t-p} + w_t + \theta_1 w_{t-1} + \cdots + \theta_q w_{t-q}$$

### 9.2. Đa thức AR và MA

- Đa thức AR: $\phi(z) = 1 - \phi_1 z - \phi_2 z^2 - \cdots - \phi_p z^p$
- Đa thức MA: $\theta(z) = 1 + \theta_1 z + \theta_2 z^2 + \cdots + \theta_q z^q$

Có thể viết gọn dưới dạng toán tử trễ $B$:

$$\phi(B) X_t = \theta(B) w_t$$

### 9.3. Điều kiện hợp lệ

1. **Tính dừng**: mọi nghiệm của $\phi(z) = 0$ nằm ngoài đường tròn đơn vị.
2. **Tính khả nghịch**: mọi nghiệm của $\theta(z) = 0$ nằm ngoài đường tròn đơn vị.
3. $\phi(z)$ và $\theta(z)$ **không có thừa số chung** (để mô hình duy nhất).

### 9.4. Bảng EACF — xác định bậc $(p, q)$

Với ARMA(p, q), cả ACF lẫn PACF đều giảm dần → khó xác định bậc trực tiếp. Ta dùng **EACF (Extended ACF)**.

Trong bảng EACF:
- `x`: hệ số tương quan mở rộng có ý nghĩa thống kê.
- `0`: không có ý nghĩa thống kê.

**Cách đọc**: tìm vị trí $(p, q)$ sao cho từ đó trở đi (sang phải và xuống dưới) đều là `0` — tạo thành **tam giác số 0**.

```r
library(TSA)
set.seed(123)
data <- sarima.sim(ar = 0.6, ma = 0.4, n = 1000)

plot(data)
adf.test(data)            # p-value < 0.05 => chuỗi dừng

par(mfrow = c(1, 2))
acf(data); pacf(data)
# Cả hai đều giảm dần => gợi ý mô hình ARMA

eacf(data)
# Tam giác 0 đỉnh tại (1, 1) => ARMA(1, 1)
```

> 📌 **EACF — Ý nghĩa & Lưu ý (bổ sung)**
>
> - **Ý nghĩa**: ACF và PACF chỉ chẩn đoán tốt với MA thuần hoặc AR thuần. Khi cả hai đều "giảm dần" (đặc trưng ARMA hỗn hợp), EACF (Tsay & Tiao, 1984) là công cụ chuẩn để dò bậc đồng thời $(p, q)$.
> - **Cách hoạt động**: lọc bỏ thành phần AR(p) trước, phần còn lại là MA(q) thuần — ACF của phần còn lại sẽ cắt cụt.
> - **Ưu điểm**: là công cụ duy nhất cho phép xác định *đồng thời* $p$ và $q$ trực quan.
> - **Nhược điểm**:
>   - Bảng EACF khá *nhiễu* trong thực tế — khó tìm tam giác 0 sạch sẽ; thường có nhiều ứng viên hợp lý.
>   - Kết quả phụ thuộc cỡ mẫu — mẫu nhỏ thì bảng không đáng tin cậy.
> - **Thực hành tốt nhất**: dùng EACF để lọc 2-3 ứng viên, sau đó so sánh AIC/BIC và chẩn đoán phần dư cho từng ứng viên — đừng tin tuyệt đối vào một vị trí.

### 9.5. Ước lượng tham số bằng `sarima()`

```r
library(astsa)

# Thử AR(2)
mod1 <- sarima(data, p = 2, d = 0, q = 0, no.constant = TRUE)
# Phương trình: x_t = 0.9361 x_{t-1} - 0.2791 x_{t-2} + w_t

# Thử ARMA(1, 1)
mod3 <- sarima(data, p = 1, d = 0, q = 1, no.constant = TRUE)
# Phương trình: x_t = 0.5609 x_{t-1} + w_t + 0.4076 w_{t-1}
```

### 9.6. Chẩn đoán mô hình

`sarima()` tự động xuất ra **4 đồ thị chẩn đoán**:

1. **Standardized residuals** — phần dư chuẩn hóa theo thời gian: kiểm tra xu hướng, phương sai thay đổi.
2. **QQ-plot phần dư** — kiểm tra tính chuẩn.
3. **ACF phần dư** — kiểm tra tự tương quan còn sót.
4. **p-value của kiểm định Ljung-Box** ở các độ trễ khác nhau.

#### Kiểm định Ljung-Box

$$Q^* = n(n+2) \sum_{k=1}^{m} \frac{\hat{\rho}_k^2}{n - k} \;\sim\; \chi^2_{m - p - q}$$

- $H_0$: **không** còn tự tương quan trong phần dư (phần dư là nhiễu trắng).
- $H_1$: còn tự tương quan.

→ Muốn mô hình tốt: **các p-value ở mọi độ trễ phải > 0.05** (chấm trên đường ngưỡng trong đồ thị).

Trong ví dụ:
- `mod1` (AR(2)): nhiều p-value < 0.05 → phần dư còn tự tương quan → mô hình **chưa đủ tốt**.
- `mod3` (ARMA(1, 1)): các p-value > 0.05, ACF phần dư trong dải xanh → phần dư **giống nhiễu trắng**.

### 9.7. So sánh mô hình bằng AIC/BIC

```r
mod2$ICs   # AIC, AICc, BIC của mô hình 2
mod3$ICs   # AIC, AICc, BIC của mô hình 3
```

AIC nhỏ hơn → mô hình tốt hơn. Trong ví dụ, `mod3$ICs` nhỏ hơn → chọn ARMA(1, 1).

### 9.8. Dự báo và khoảng tin cậy

```r
dubao <- sarima.for(data, n.ahead = 10, p = 1, d = 0, q = 1)

dubao$pred                                      # giá trị dự báo
dubao$se                                        # sai số chuẩn của dự báo

canduoi <- dubao$pred - 1.96 * dubao$se         # cận dưới 95%
cantren <- dubao$pred + 1.96 * dubao$se         # cận trên 95%
```

Khoảng tin cậy 95% cho dự báo:

$$\widehat{X}_{n+h} \pm 1.96 \cdot \mathrm{SE}(\widehat{X}_{n+h})$$

> 📌 **ARMA(p,q) & Kiểm định Ljung-Box — Bổ sung**
>
> **ARMA(p, q)**:
> - **Ý nghĩa**: kết hợp sức mạnh của AR (mô tả phụ thuộc dài hạn) và MA (mô tả các cú sốc gần đây) trong cùng một mô hình.
> - **Ưu điểm**:
>   - Tiết kiệm tham số: thường ARMA bậc thấp đã đủ thay thế cho AR hoặc MA bậc cao tương đương.
>   - Linh hoạt biểu diễn nhiều dạng phụ thuộc.
> - **Nhược điểm**:
>   - Yêu cầu chuỗi dừng (cần ARIMA cho chuỗi không dừng).
>   - Mô hình **tuyến tính** — không bắt được sốc đột ngột, thay đổi chế độ, hiệu ứng phi tuyến.
>   - Ràng buộc tham số phức tạp (cả dừng + khả nghịch + không có thừa số chung).
>   - Việc xác định bậc khá *chủ quan* — hai chuyên gia có thể chọn khác nhau.
> - **Khi nào nên dùng**: chuỗi dừng, không có mùa vụ, không có sốc đột ngột, dự báo ngắn-trung hạn.
> - **Khi nào KHÔNG nên dùng**: chuỗi tài chính có biến động cụm (cần GARCH); chuỗi có chế độ chuyển (cần TAR/MS-AR); chuỗi có mùa vụ (cần SARIMA).
>
> **Kiểm định Ljung-Box** (chi tiết):
> - **Thống kê**: $Q^* = n(n+2) \sum_{k=1}^{m} \dfrac{\hat{\rho}_k^2}{n-k}$ — kiểm tra **tổng thể** xem ACF của phần dư tại các độ trễ $1, \ldots, m$ có khác 0 đáng kể không.
> - **Bậc tự do**: với phần dư của ARIMA(p,d,q), bậc tự do = $m - p - q$ (không phải $m$).
> - **Quy tắc chọn $m$**: $m = \min(10, n/5)$ với chuỗi không mùa; $m = \min(2s, n/5)$ với chuỗi có chu kỳ $s$.
> - **Mong muốn**: *không* bác bỏ $H_0$ → mọi $p\text{-value} > 0.05$ → phần dư là nhiễu trắng.
> - **Hạn chế**: Là **portmanteau test** — chỉ phát hiện tự tương quan trung bình; có thể bỏ sót tự tương quan tại một độ trễ riêng lẻ.

---

# Phần IV — Mô hình ARIMA & SARIMA

## Chương 10. Mô hình ARIMA(p, d, q)

> Nguồn: [Buoi25_3.R](notebooks/Buoi25_3.R)

### 10.1. Quy trình xây dựng mô hình (Box-Jenkins)

1. **Vẽ biểu đồ** chuỗi gốc — quan sát xu hướng, mùa vụ, phương sai.
2. **Biến đổi** dữ liệu nếu cần:
   - Lấy `log()` để ổn định phương sai.
   - Lấy `diff()` để loại bỏ xu hướng.
   - Biến đổi lũy thừa Box-Cox.
3. **Xác định bậc** $(p, d, q)$ qua ACF, PACF, EACF.
4. **Ước lượng tham số**.
5. **Chẩn đoán & lựa chọn** mô hình (phần dư + AIC/BIC).
6. **Dự báo**.

### 10.2. Định nghĩa ARIMA(p, d, q)

$$\phi(B) (1 - B)^d X_t = \theta(B) w_t$$

trong đó $d$ là số lần sai phân để chuỗi trở thành dừng.

### 10.3. Ví dụ mô phỏng ARIMA(1, 1, 1)

```r
library(astsa)
library(tseries)

set.seed(123)
data <- sarima.sim(n = 1000, ar = 0.5, d = 1, ma = 0.2)
plot(data)

adf.test(data)
# p-value = 0.3688 > 0.05 => chuỗi không dừng

adf.test(diff(data))
# p-value = 0.01 < 0.05 => chuỗi sau sai phân là dừng

par(mfrow = c(1, 2))
acf(diff(data))     # giảm dần
pacf(diff(data))    # gợi ý AR(2)

library(TSA)
eacf(diff(data))    # gợi ý ARMA(1, 1) hoặc MA(3)
```

Fit và so sánh ba ứng viên:

```r
mod1 <- sarima(diff(data), p = 2, d = 0, q = 0, no.constant = TRUE)  # AR(2)
mod2 <- sarima(diff(data), p = 1, d = 0, q = 1, no.constant = TRUE)  # ARMA(1, 1)
mod3 <- sarima(diff(data), p = 0, d = 0, q = 3, no.constant = TRUE)  # MA(3)

# Phần dư cả 3 mô hình đều giống nhiễu trắng
mod1$ICs; mod2$ICs; mod3$ICs
# Mô hình có AIC nhỏ nhất là tốt nhất
```

### 10.4. Ví dụ thực tế với dữ liệu `deere3`

```r
library(TSA)
data("deere3")
plot(deere3)

adf.test(deere3)
# p-value = 0.01 < 0.05 => chuỗi dừng

acf(deere3)    # cắt ở độ trễ 2
pacf(deere3)   # cắt ở độ trễ 1

ar1 <- sarima(deere3, p = 1, d = 0, q = 0, no.constant = TRUE)  # AR(1)
ar2 <- sarima(deere3, p = 2, d = 0, q = 0, no.constant = TRUE)  # AR(2)
```

Cả hai mô hình:
- ACF của phần dư không có giá trị nào lớn đáng kể.
- $p\text{-value}$ kiểm định Ljung-Box đều > 0.05 → phần dư không còn tự tương quan.
- Các điểm trên QQ-plot khá sát đường thẳng chuẩn → phần dư có tính chuẩn.

→ Cả hai mô hình đều có phần dư giống nhiễu trắng.

```r
ar1$ICs   # nhỏ hơn => AR(1) tốt hơn theo AIC
ar2$ICs

sarima.for(deere3, n.ahead = 10, p = 1, d = 0, q = 0)
```

### 10.5. Tự động chọn bậc với `auto.arima()`

```r
install.packages("forecast")
library(forecast)
auto.arima(deere3)
```

Hàm `auto.arima()` của package **forecast** quét nhiều tổ hợp $(p, d, q)$ và chọn mô hình theo AICc.

> Lưu ý: `auto.arima()` là công cụ hỗ trợ — kết quả nên được kiểm chứng lại bằng chẩn đoán phần dư và so sánh với mô hình thủ công.

> 📌 **ARIMA — Khái niệm, Ý nghĩa, Ưu/Nhược (bổ sung)**
>
> - **Khái niệm**: ARIMA = AR + Integrated (sai phân) + MA. Tham số $d$ cho biết số lần sai phân để chuỗi trở thành dừng.
> - **Ý nghĩa**: là *bộ khung tổng quát* của các mô hình tuyến tính cho chuỗi thời gian không dừng — mở rộng ARMA cho dữ liệu thực tế.
> - **Ưu điểm**:
>   - Linh hoạt: mô tả được xu hướng, dao động, phụ thuộc dài hạn trong một khung thống nhất.
>   - Dựa trên nền tảng thống kê chặt chẽ → cho khoảng tin cậy hợp lệ.
>   - **Hiệu quả với chuỗi vừa và nhỏ** — chỉ cần dữ liệu chuỗi đó, không cần biến ngoại lai.
>   - Hoạt động tốt cho **dự báo ngắn hạn**.
> - **Nhược điểm**:
>   - **Tuyến tính** — bất lực với sốc đột ngột, regime change, hiệu ứng phi tuyến.
>   - **Khó dự đoán điểm uốn** (turning points).
>   - Việc chọn $(p, d, q)$ mang tính *chủ quan*, đòi hỏi kinh nghiệm.
>   - **Kém với dự báo dài hạn** — sai số tích lũy.
>   - Nhạy cảm với outliers và missing values.
>   - Không tích hợp được biến giải thích bên ngoài (cần ARIMAX/Regression with ARIMA errors).
>   - **Không xử lý được mùa vụ** trực tiếp — cần SARIMA.
> - **Khi nào nên dùng**: chuỗi đơn biến, không có mùa vụ mạnh, dữ liệu sạch, dự báo ngắn-trung hạn (vài bước → vài chục bước).
> - **Khi nào KHÔNG nên dùng**:
>   - Chuỗi quá ngắn (< 50 quan sát) hoặc quá dài có thay đổi cấu trúc.
>   - Chuỗi có mùa vụ → dùng SARIMA / TBATS.
>   - Chuỗi có biến độc lập quan trọng → dùng ARIMAX hoặc mô hình hồi quy.
>   - Chuỗi có biến động cụm (volatility clustering) → dùng GARCH.
>   - Dự báo nhiều bước về tương lai (cần mô hình state-space hoặc deep learning).

> 📌 **`auto.arima()` — Khi nào tin, khi nào kiểm tra lại (bổ sung)**
>
> - **Cơ chế**: thuật toán Hyndman-Khandakar — kết hợp kiểm định unit root + tối thiểu AICc + MLE. Mặc định dùng *stepwise search* (nhanh nhưng có thể bỏ sót mô hình tốt).
> - **Ưu điểm**: tiết kiệm thời gian, không cần am hiểu sâu Box-Jenkins, phù hợp khi cần xử lý hàng loạt nhiều chuỗi.
> - **Nhược điểm**:
>   - Mô hình AICc thấp nhất *trên tập huấn luyện* không đảm bảo tốt nhất *trên tập kiểm tra*.
>   - Tìm kiếm stepwise có thể bỏ qua mô hình tối ưu.
>   - Không thay thế được kinh nghiệm chuyên gia khi chuỗi có đặc thù bất thường.
> - **Khuyến nghị thực hành**:
>   - Khi phân tích **một chuỗi đơn lẻ**: bật `stepwise = FALSE, approximation = FALSE` để quét đầy đủ.
>   - Luôn **chẩn đoán phần dư** của mô hình `auto.arima` chọn ra — không tin mù quáng.
>   - So sánh với mô hình thủ công dựa trên ACF/PACF/EACF.

---

## Chương 11. Phân rã chuỗi thời gian

> Nguồn: [Buoi1_4.R](notebooks/Buoi1_4.R)

### 11.1. Mô hình cộng tính

$$X_t = T_t + S_t + R_t$$

- $T_t$: thành phần xu hướng.
- $S_t$: thành phần mùa vụ.
- $R_t$: phần dư (random / remainder).

```r
decompose(co2)
plot(decompose(co2))
```

### 11.2. Mô hình nhân tính

$$X_t = T_t \cdot S_t \cdot R_t$$

Phù hợp khi biên độ dao động mùa vụ tăng theo mức của chuỗi (ví dụ chuỗi `airpass` — số hành khách hàng không).

```r
library(TSA)
plot(airpass)
decompose(airpass, type = "multiplicative")
plot(decompose(airpass, type = "multiplicative"))
```

### 11.3. Trung bình trượt (moving average)

```r
library(forecast)
co2_ma <- ma(co2, order = 12)
plot(co2)
lines(co2_ma, col = "red")
```

`ma(x, order = k)` tính trung bình trượt với cửa sổ kích thước $k$ — dùng để **làm mượt** chuỗi và làm nổi bật xu hướng.

> 📌 **Phân rã chuỗi thời gian — Bổ sung**
>
> - **Ý nghĩa**: tách chuỗi thành các thành phần ý nghĩa kinh tế/vật lý — phục vụ phân tích mô tả, hiệu chỉnh mùa vụ (seasonally adjusted series), và làm tiền đề cho dự báo.
> - **Mô hình cộng tính** ($X_t = T_t + S_t + R_t$): dùng khi **biên độ mùa vụ KHÔNG đổi** theo mức của chuỗi.
> - **Mô hình nhân tính** ($X_t = T_t \cdot S_t \cdot R_t$): dùng khi **biên độ mùa vụ TĂNG/GIẢM theo mức** của chuỗi (rất phổ biến với chuỗi tài chính, hành khách hàng không như `airpass`).
>   - *Mẹo*: lấy $\log$ chuỗi nhân tính sẽ chuyển thành chuỗi cộng tính.
> - **Ưu điểm chung của `decompose()` cổ điển**: đơn giản, dễ giải thích, tính nhanh.
> - **Nhược điểm**:
>   - **Mất dữ liệu ở 2 đầu chuỗi** (trung bình trượt cần $k$ điểm xung quanh).
>   - Giả định mùa vụ *lặp lại y hệt* mỗi chu kỳ → kém với mùa vụ biến đổi dần.
>   - **Nhạy cảm với outlier**.
>   - Không cung cấp khoảng tin cậy.
> - **Thay thế hiện đại**: STL (Seasonal-Trend decomposition with Loess) cho phép mùa vụ thay đổi theo thời gian, robust với outlier; X-13ARIMA-SEATS cho dữ liệu kinh tế chính thức.

---

## Chương 12. Mô hình SARIMA(p, d, q)(P, D, Q)_s

### 12.1. Định nghĩa

$$\Phi(B^s)\, \phi(B)\, (1 - B)^d (1 - B^s)^D X_t = \Theta(B^s)\, \theta(B)\, w_t$$

trong đó:
- $(p, d, q)$: bậc AR, sai phân, MA cho thành phần **không mùa**.
- $(P, D, Q)$: bậc AR, sai phân, MA cho thành phần **mùa vụ** với chu kỳ $s$.
- $s$: chu kỳ mùa (ví dụ $s = 12$ với dữ liệu tháng có chu kỳ năm).

### 12.2. Ví dụ với dữ liệu `co2`

```r
plot(co2)                 # nhận thấy có tính mùa khá rõ rệt

library(tseries)
adf.test(co2)
# p-value < 0.05 => chuỗi dừng (theo nghĩa không mùa)
# => chỉ cần sai phân theo mùa

diff_co2 <- diff(co2, lag = 12)
plot(diff_co2)
```

### 12.3. Phân tích ACF/PACF với `lag.max`

```r
par(mfrow = c(1, 2))
acf(diff(co2, lag = 12), lag.max = 36)
pacf(diff(co2, lag = 12), lag.max = 36)
```

Quan sát:
- ACF giảm dần.
- PACF có giá trị đáng kể tại các độ trễ 1, 2, 12 và 24.

→ Dự đoán bậc:
- $(p, d, q) = (2, 0, 0)$
- $(P, D, Q) = (2, 1, 0)$, $s = 12$.
- Mô hình: **ARIMA(2, 0, 0) × (2, 1, 0)$_{12}$**.

### 12.4. Fit mô hình SARIMA

```r
library(astsa)
mod1 <- sarima(co2,
               p = 2, d = 0, q = 0,
               P = 2, D = 1, Q = 0,
               S = 12,
               no.constant = TRUE)
```

### 12.5. So sánh với `auto.arima()`

```r
library(forecast)
auto.arima(co2)
# Đề xuất: ARIMA(1, 0, 1)(0, 1, 1)_12

mod2 <- sarima(co2,
               p = 1, d = 0, q = 1,
               P = 0, D = 1, Q = 1, S = 12)

mod1$ICs
mod2$ICs
# Mô hình có AIC nhỏ hơn được ưu tiên
```

### 12.6. Dự báo dài hạn

```r
par(mfrow = c(1, 1))
sarima.for(co2, n.ahead = 36, p = 1, d = 0, q = 1,
           P = 0, D = 1, Q = 1, S = 12)
```

Hàm `sarima.for()` tự động vẽ chuỗi gốc + chuỗi dự báo + dải tin cậy 95%.

### 12.7. Kết hợp xu hướng + mùa bằng hồi quy

Thay vì SARIMA, một cách tiếp cận khác là kết hợp `season()` (mùa vụ) và `time()` (xu hướng tuyến tính):

```r
mod_reg <- lm(co2 ~ season(co2) + time(co2))
summary(mod_reg)
acf(residuals(mod_reg))
```

→ Quan sát ACF của phần dư để đánh giá mức độ phù hợp của mô hình hồi quy này.

> 📌 **SARIMA — Khái niệm, Ý nghĩa, Ưu/Nhược (bổ sung)**
>
> - **Khái niệm**: mở rộng ARIMA bằng cách thêm các thành phần AR, MA, sai phân **theo bội số chu kỳ mùa $s$**. Sai phân mùa $(1 - B^{12}) X_t = X_t - X_{t-12}$ loại bỏ thành phần mùa.
> - **Ý nghĩa**: là *vũ khí chuẩn* cho chuỗi có cả xu hướng + mùa vụ — kết hợp được phụ thuộc trong cùng mùa và phụ thuộc giữa các mùa kề nhau.
> - **Ưu điểm**:
>   - Bắt được mùa vụ **tự nhiên** mà ARIMA thuần không xử lý được.
>   - Tiết kiệm tham số hơn so với việc dùng ARIMA bậc cao ($p \geq s$) cố mô phỏng mùa.
>   - Có nền tảng lý thuyết chặt → khoảng tin cậy hợp lệ.
> - **Nhược điểm**:
>   - **Phức tạp**: 7 tham số $(p, d, q, P, D, Q, s)$ → khó chọn, cần kinh nghiệm.
>   - **Yêu cầu dữ liệu dài**: cần ít nhất 3-4 chu kỳ mùa đầy đủ để ước lượng tin cậy.
>   - **Dễ overfit** với chuỗi ngắn hoặc nhiều tham số.
>   - **Không xử lý chu kỳ mùa rất dài** ($s > 200$): gặp vấn đề bộ nhớ và bất ổn số học → nên dùng hồi quy điều hòa hoặc TBATS.
>   - **Một mùa vụ duy nhất**: không phù hợp khi chuỗi có nhiều chu kỳ mùa cùng lúc (ví dụ điện năng có chu kỳ ngày + tuần + năm).
>   - **Mùa vụ cố định**: SARIMA giả định pattern mùa không thay đổi.
> - **Khi nào nên dùng**: chuỗi có mùa vụ rõ rệt, chu kỳ ngắn ($s \leq 24$), dữ liệu đủ dài (≥ 3 chu kỳ), pattern mùa ổn định.
> - **Khi nào KHÔNG nên dùng**:
>   - Chu kỳ mùa rất dài → dùng *Dynamic Harmonic Regression* hoặc TBATS.
>   - Nhiều mùa vụ chồng chất → dùng TBATS, Prophet.
>   - Mùa vụ biến đổi nhiều → dùng STL + ETS hoặc Prophet.
>   - Chuỗi quá ngắn (< 2-3 chu kỳ) → dùng phương pháp đơn giản hơn (Holt-Winters, decomposition + naive).

---

# Phần V — Phụ lục

## A. Tổng hợp các package R sử dụng

| Package | Hàm chính | Chức năng |
|---------|-----------|-----------|
| `tseries` | `adf.test()`, `kpss.test()` | Kiểm định tính dừng |
| `TSA` | `season()`, `harmonic()`, `zlag()`, `eacf()`, `runs()` | Hồi quy xu hướng, EACF, run test, datasets |
| `astsa` | `sarima()`, `sarima.sim()`, `sarima.for()` | Fit ARIMA/SARIMA, mô phỏng, dự báo |
| `forecast` | `auto.arima()`, `ma()` | Tự động chọn mô hình, trung bình trượt |
| `stats` (base) | `lm()`, `decompose()`, `acf()`, `pacf()`, `diff()`, `predict()`, `shapiro.test()`, `qqnorm()`, `qqline()` | Hồi quy, phân rã, ACF/PACF, biến đổi cơ bản |

## B. Tổng hợp các kiểm định thống kê

| Kiểm định | Mục đích | $H_0$ | Bác bỏ $H_0$ khi | Hàm R |
|-----------|----------|-------|---------|-------|
| **ADF** | Tính dừng | Chuỗi không dừng | $p < 0.05$ ⇒ **dừng** | `adf.test()` (`tseries`) |
| **KPSS** | Tính dừng | Chuỗi dừng | $p < 0.05$ ⇒ **không dừng** | `kpss.test()` (`tseries`) |
| **Shapiro-Wilk** | Tính chuẩn | Phân phối chuẩn | $p < 0.05$ ⇒ **không chuẩn** | `shapiro.test()` (base) |
| **Run test** | Tính độc lập / ngẫu nhiên | Chuỗi ngẫu nhiên | $p < 0.05$ ⇒ **có phụ thuộc** | `runs()` (`TSA`) |
| **Ljung-Box** | Tự tương quan tổng quát | Không tự tương quan | $p < 0.05$ ⇒ **còn tự tương quan** | tự động trong `sarima()` (`astsa`) |

## C. Bảng nhận diện ACF/PACF cho MA / AR / ARMA

| Mô hình | ACF | PACF | Công cụ chính |
|---------|-----|------|---------------|
| **MA(q)** | **Cắt cụt tại độ trễ $q$** | Giảm dần | ACF |
| **AR(p)** | Giảm dần | **Cắt cụt tại độ trễ $p$** | PACF |
| **ARMA(p, q)** | Giảm dần | Giảm dần | **EACF** (tam giác số 0) |
| **ARIMA(p, d, q)** | Giảm rất chậm khi chưa sai phân | — | Sai phân + ACF/PACF + EACF |
| **SARIMA** | Có đỉnh tại độ trễ bội số chu kỳ $s$ | Có đỉnh tại độ trễ bội số $s$ | Sai phân mùa + ACF/PACF với `lag.max ≥ 2s` |

## D. Sơ đồ quy trình Box-Jenkins

> 📌 **Phương pháp luận Box-Jenkins — Bổ sung**
>
> - **Khái niệm**: do George Box & Gwilym Jenkins đề xuất (1970), là quy trình lặp gồm 3 bước cốt lõi: **Identification → Estimation → Diagnostic Checking**, sau đó là *Forecasting*.
> - **Ý nghĩa**: cung cấp một **khung tư duy có hệ thống** để xây dựng mô hình ARIMA — từ trực giác đồ thị đến kiểm tra thống kê chặt chẽ.
> - **Ưu điểm**:
>   - Rõ ràng, có thể lặp lại; xử lý được cả xu hướng, mùa vụ và phụ thuộc.
>   - Phù hợp với dự báo ngắn hạn và phân tích kinh tế-tài chính.
> - **Nhược điểm**:
>   - **Tốn thời gian**, đòi hỏi kinh nghiệm.
>   - **Chủ quan ở bước Identification**: hai người dùng có thể chọn hai mô hình khác nhau cho cùng dữ liệu.
>   - **Cần dữ liệu lịch sử nhiều** (thường ≥ 50 quan sát).
>   - Kém với chuỗi rất nhiễu hoặc có pattern bất thường.
> - **Khi nào nên dùng**: chuỗi đơn biến có cấu trúc thống kê ổn định, độ dài vừa-lớn, không có sốc cấu trúc lớn.
> - **Khi nào KHÔNG nên dùng**: chuỗi quá ngắn; chuỗi với biến giải thích quan trọng (dùng hồi quy + ARIMA errors); chuỗi phi tuyến / có chế độ chuyển (dùng TAR, MS-AR, ML).

```
┌──────────────────────────────────────────────────────────┐
│  1. THĂM DÒ DỮ LIỆU                                      │
│     - plot(x)                                            │
│     - Quan sát xu hướng, mùa vụ, phương sai              │
└──────────────────────┬───────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────┐
│  2. KIỂM TRA TÍNH DỪNG                                   │
│     - adf.test(x)   (H0: không dừng)                     │
│     - kpss.test(x)  (H0: dừng)                           │
└──────────────────────┬───────────────────────────────────┘
                       │
              ┌────────┴────────┐
              ▼                 ▼
         CHƯA DỪNG           ĐÃ DỪNG
              │                 │
              ▼                 │
┌──────────────────────────┐    │
│  3. BIẾN ĐỔI             │    │
│     - log(x)             │    │
│     - diff(x), diff(x,12)│    │
│     - Box-Cox            │    │
└──────────┬───────────────┘    │
           │                    │
           └─────────┬──────────┘
                     ▼
┌──────────────────────────────────────────────────────────┐
│  4. XÁC ĐỊNH BẬC                                         │
│     - acf(x), pacf(x)                                    │
│     - eacf(x)                                            │
│     - auto.arima(x)                                      │
└──────────────────────┬───────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────┐
│  5. ƯỚC LƯỢNG THAM SỐ                                    │
│     - sarima(x, p, d, q[, P, D, Q, S])                   │
└──────────────────────┬───────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────┐
│  6. CHẨN ĐOÁN PHẦN DƯ                                    │
│     - Standardized residuals plot                        │
│     - QQ-plot + Shapiro-Wilk  (tính chuẩn)               │
│     - ACF của phần dư                                    │
│     - Ljung-Box (tất cả p-value > 0.05?)                 │
└──────────────────────┬───────────────────────────────────┘
                       │
              ┌────────┴────────┐
              ▼                 ▼
         CHƯA TỐT             TỐT
              │                 │
       quay lại bước 4          ▼
                       ┌──────────────────────────────────┐
                       │  7. SO SÁNH MÔ HÌNH (AIC / BIC)  │
                       │     - mod$ICs                    │
                       └──────────────┬───────────────────┘
                                      │
                                      ▼
                       ┌──────────────────────────────────┐
                       │  8. DỰ BÁO                       │
                       │     - sarima.for()               │
                       │     - Khoảng tin cậy 95%         │
                       └──────────────────────────────────┘
```

## E. Danh mục datasets sử dụng

| Dataset | Package | Mô tả ngắn | Sử dụng tại chương |
|---------|---------|------------|--------------------|
| `oil.price` | `TSA` | Giá dầu | Ch. 1, 2 |
| `rwalk` | `TSA` | Chuỗi bước ngẫu nhiên (random walk) | Ch. 4 |
| `tempdub` | `TSA` | Nhiệt độ trung bình hàng tháng tại Dubuque | Ch. 5, 6 |
| `ma1.2.s` | `TSA` | Chuỗi mô phỏng MA(1) | Ch. 7 |
| `ma2.s` | `TSA` | Chuỗi mô phỏng MA(2) | Ch. 7 |
| `ar1.2.s` | `TSA` | Chuỗi mô phỏng AR(1) | Ch. 8 |
| `ar2.s` | `TSA` | Chuỗi mô phỏng AR(2) | Ch. 8 |
| `soi` | `astsa` | Southern Oscillation Index | Ch. 8 |
| `deere3` | `TSA` | Sai số gia công Deere | Ch. 10 |
| `co2` | `datasets` (base) | Nồng độ CO₂ tại Mauna Loa | Ch. 11, 12 |
| `airpass` | `TSA` | Số hành khách hàng không quốc tế | Ch. 11 |

---

> **Cách sử dụng báo cáo này**:
> 1. Đọc theo thứ tự để nắm hệ thống lý thuyết.
> 2. Mở song song file `.R` tương ứng (đường link đầu mỗi chương) để chạy lại code.
> 3. Đối chiếu với phụ lục B & C khi cần tra cứu nhanh kiểm định / đặc trưng mô hình.

---

## Tài liệu tham khảo bổ sung (từ web)

**Tính dừng — ADF & KPSS**
- [Stationarity and detrending (ADF/KPSS) — statsmodels](https://www.statsmodels.org/stable/examples/notebooks/generated/stationarity_detrending_adf_kpss.html)
- [Comparison study of ADF vs KPSS — Tannya Sharma, Medium](https://medium.com/@tannyasharma21/comparision-study-of-adf-vs-kpss-test-c9d8dec4f62a)
- [KPSS Test for Stationarity — Machine Learning Plus](https://www.machinelearningplus.com/time-series/kpss-test-for-stationarity/)

**Biến đổi chuỗi — sai phân, log**
- [Differencing and Log Transformation — Finance Train](https://financetrain.com/differencing-and-log-transformation)
- [Stationarity and differencing — Hyndman, FPP2 §8.1](https://otexts.com/fpp2/stationarity.html)
- [Uses of the logarithm transformation — Robert Nau, Duke](https://people.duke.edu/~rnau/411log.htm)

**Xu hướng tất định vs ngẫu nhiên**
- [Stochastic and deterministic trends — Hyndman, FPP3 §10.4](https://otexts.com/fpp3/stochastic-and-deterministic-trends.html)
- [Trend-Stationary vs Difference-Stationary — MathWorks](https://www.mathworks.com/help/econ/trend-stationary-vs-difference-stationary.html)

**Hồi quy điều hòa (Fourier / cosine)**
- [Dynamic harmonic regression — Hyndman, FPP3 §10.5](https://otexts.com/fpp3/dhr.html)
- [What is Harmonic Regression for Time Series — Towards Data Science](https://towardsdatascience.com/take-your-forecasting-to-the-next-level-with-harmonic-regression-5a8515f63295/)

**AR, MA, ARMA**
- [Building blocks of time series modelling: AR, MA & ARMA — Medium](https://medium.com/@div09/building-blocks-of-time-series-modelling-ar-ma-arma-explained-b7889733da92)
- [Rules for identifying ARIMA models — Robert Nau, Duke](https://people.duke.edu/~rnau/arimrule.htm)
- [Autoregressive moving-average model — Wikipedia](https://en.wikipedia.org/wiki/Autoregressive_moving-average_model)

**EACF**
- [eacf: Compute the sample ESACF — TSA package R](https://rdrr.io/cran/TSA/man/eacf.html)
- [Extended Autocorrelation Function — freshrimpsushi](https://freshrimpsushi.github.io/en/posts/1213/)

**ARIMA**
- [Understanding the Limitations of ARIMA Forecasting — Towards Data Science](https://towardsdatascience.com/understanding-the-limitations-of-arima-forecasting-899f8d8e5cf3/)
- [Understanding ARIMA models for machine learning — Capital One](https://www.capitalone.com/tech/machine-learning/understanding-arima-models/)
- [ARIMA — Wikipedia](https://en.wikipedia.org/wiki/Autoregressive_integrated_moving_average)

**SARIMA**
- [Seasonal ARIMA models — Hyndman, FPP3 §9.9](https://otexts.com/fpp3/seasonal-arima.html)
- [SARIMA: Complete Guide — Michael Brenndoerfer](https://mbrenndoerfer.com/writing/sarima-seasonal-time-series-forecasting)
- [ARIMA vs SARIMA — GeeksforGeeks](https://www.geeksforgeeks.org/machine-learning/arima-vs-sarima-model/)

**Box-Jenkins**
- [Box–Jenkins method — Wikipedia](https://en.wikipedia.org/wiki/Box%E2%80%93Jenkins_method)
- [Box-Jenkins Methodology — Columbia Mailman School of Public Health](https://www.publichealth.columbia.edu/research/population-health-methods/box-jenkins-methodology)
- [A Gentle Introduction to the Box-Jenkins Method — Machine Learning Mastery](https://machinelearningmastery.com/gentle-introduction-box-jenkins-method-time-series-forecasting/)

**Kiểm định Ljung-Box**
- [Ljung–Box test — Wikipedia](https://en.wikipedia.org/wiki/Ljung%E2%80%93Box_test)
- [Thoughts on the Ljung-Box test — Rob J Hyndman](https://robjhyndman.com/hyndsight/ljung-box-test/)
- [Ljung-Box Test: Definition + Example — Statology](https://www.statology.org/ljung-box-test/)

**auto.arima()**
- [auto.arima reference — forecast package, Hyndman](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
- [ARIMA modelling in R — Hyndman, FPP2 §8.7](https://otexts.com/fpp2/arima-r.html)

**Phân rã chuỗi thời gian**
- [Classical decomposition — Hyndman, FPP3 §3.4](https://otexts.com/fpp3/classical-decomposition.html)
- [Decomposition of time series — Wikipedia](https://en.wikipedia.org/wiki/Decomposition_of_time_series)
- [Different Types of Time Series Decomposition — Towards Data Science](https://towardsdatascience.com/different-types-of-time-series-decomposition-396c09f92693/)
