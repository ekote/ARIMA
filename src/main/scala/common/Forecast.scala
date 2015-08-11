package common

import breeze.linalg.{ DenseVector, DenseMatrix, reverse, sum }
import breeze.numerics._
import breeze.stats.{ mean }
import common.ARIMAUtils.{ diff, diffinv, ma_inf, smooth_ma }
import common.Innovations._

object Forecast {
  private def forecast_arma(y: DenseVector[Double], model: (DenseVector[Double], DenseVector[Double], Double), h: Int, demean: Boolean): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = y.length
    var mu = mean(y)
    if (!demean) mu = 0.0
    var x = y - mu
    val phi = model._1
    val theta = model._2
    var dx = innovations(x = x, model = model)
    dx = DenseVector.vertcat(dx, DenseVector.zeros[Double](h))
    x = DenseVector.vertcat(x, DenseVector.zeros[Double](h))
    val p = phi.length
    val q = theta.length
    for (t <- n to n + h - 1) {
      val A = phi dot reverse(x(t - p to t - 1))
      val B = theta dot reverse(dx(t - q to t - 1))
      x(t) = A + B
    }
    val pred = x(n to n + h - 1) + mu
    (pred, phi, null, null)
  }

  private def forecast_diff(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), d: Int, h: Int, k: Int, demean: Boolean): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = x.length
    var lag = 1
    if (xv(k - 1) == "diff") lag = xv(k).toInt
    val differences = d
    val y = diff(x, lag, differences)
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, d = d, h = h, k = k + 2, demean = demean)
    pred = diffinv(pred, lag, differences, x(n - lag * differences to n - 1))
    pred = pred(lag to lag + h - 1)
    if (phi == null) sys.error("diff before log")
    phi = DenseVector[Double](1.0 +: (-phi).toArray)
    phi = DenseVector.vertcat(phi, DenseVector.fill[Double](lag, 0.0)) - DenseVector.vertcat(DenseVector.fill[Double](lag, 0.0), phi)
    phi = -phi(1 until phi.length)
    (pred, phi, l, u)
  }

  private def forecast_log(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), d: Int, h: Int, k: Int, demean: Boolean): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    var (pred, phi, l, u) = forecast_transform(x = log(x), xv = xv, model = model, d = d, h = h, k = k + 1, demean = demean)
    if (phi != null) {
      val theta = model._2
      val sigma2 = model._3
      val psi = ma_inf(phi = phi, theta = theta, n = h)
      def g(j: Int) = sum(pow(psi(0 to j - 1), 2))
      val se = DenseVector.zeros[Double](h)
      for (i <- 1 to h) se(i - 1) = sqrt(sigma2 * g(i))
      l = pred - (se * 1.96)
      u = pred + (se * 1.96)
    }
    pred = exp(pred)
    l = exp(l)
    u = exp(u)
    (pred, null, l, u)
  }

  private def forecast_season(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), d: Int, h: Int, k: Int, demean: Boolean): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = x.length
    val d = xv(k).toInt
    val q = d / 2
    def F1(t: Int): Double = x(t - q - 1) / 2 + sum(x(t - q to t + q - 2)) + x(t + q - 1) / 2
    var m = DenseVector.zeros[Double](n - 2 * q)
    for (i <- q + 1 to n - q) m(i - q - 1) = F1(i) / d
    m = DenseVector[Double](Array.fill[Double](q)(0) ++ m.toArray ++ Array.fill[Double](q)(0))
    if (d != 2 * q) m = smooth_ma(x, q)
    val dx = x - m
    def F2(k: Int): Double = mean(dx(k + q - 1 to n - q - 1 by d))
    var w = DenseVector.zeros[Double](d)
    for (i <- 1 to d) w(i - 1) = F2(i)
    w = w - mean(w)
    var s = w
    for (i <- 1 until n + q + h) s = DenseVector.vertcat(s, w)
    s = s(q to q + n + h - 1)
    val y = x - s(0 until n)
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, d = d, h = h, k = k + 2, demean = demean)
    pred = pred + s(n to n + h - 1)
    if (phi == null) {
      l = l + s(n to n + h - 1)
      u = u + s(n to n + h - 1)
    }
    (pred, phi, l, u)
  }

  private def forecast_trend(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), d: Int, h: Int, k: Int = 1, demean: Boolean): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = x.length
    val p = xv(k).toInt
    val X = DenseMatrix.zeros[Double](n + h, p + 1)
    for (j <- 0 to p) X(::, j) := pow(DenseVector.rangeD(1.0, n + h + 1), j)
    val b = X(0 until n, ::) \ x
    val xhat = X * b
    val y = x - xhat(0 until n)
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, d = d, h = h, k = k + 2, demean = demean)
    pred = pred + xhat(n until n + h)
    if (phi == null) {
      l = l + xhat(n until n + h)
      u = u + xhat(n until n + h)
    }
    (pred, phi, l, u)
  }

  private def forecast_transform(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), d: Int, h: Int, k: Int = 1, demean: Boolean): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    if (k > xv.length) return forecast_arma(y = x, model = model, h = h, demean = demean)
    if (xv(k - 1) == "diff") return forecast_diff(x = x, xv = xv, model = model, d = d, h = h, k = k, demean = demean)
    if (xv(k - 1) == "log") return forecast_log(x = x, xv = xv, model = model, d = d, h = h, k = k, demean = demean)
    if (xv(k - 1) == "season") return forecast_season(x = x, xv = xv, model = model, d = d, h = h, k = k, demean = demean)
    if (xv(k - 1) == "trend") return forecast_trend(x = x, xv = xv, model = model, d = d, h = h, k = k, demean = demean)
    else sys.error("x vector transformation is invalid")
  }

  def forecast(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), d: Int, h: Int, k: Int = 1, demean: Boolean, enableBound: Boolean = false, leastBound: Int, upperBound: Int): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    var (pred, phi, l, u) = forecast_transform(x = x, xv = xv, model = model, d = d, h = h, k = 1, demean = demean)
    val se = DenseVector.zeros[Double](h)
    if (enableBound) pred = pred.map(v => if (v < leastBound) leastBound else if (v > upperBound) upperBound else v)
    if (phi != null) {
      val theta = model._2
      val psi = ma_inf(phi = phi, theta = theta, n = h)
      val sigma2 = model._3
      def g(j: Int): Double = sum(pow(psi(0 until j), 2))
      for (i <- 1 to h) se(i - 1) = sqrt(sigma2 * g(i))
      l = pred - (se * 1.96)
      u = pred + (se * 1.96)
    } else return (pred, null, l, u)

    (pred, se, l, u)
  }
}