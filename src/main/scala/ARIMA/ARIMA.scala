
/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

// package TimeSeriesForecast.ARIMA

package ARIMA

import breeze.linalg.{DenseVector, DenseMatrix, max, inv, diag, sum, all, reverse}
import breeze.numerics.{abs, log, NaN, pow, sqrt, exp}
import breeze.optimize.{ApproximateGradientFunction, LBFGS, DiffFunction}
import breeze.stats.mean
import common.ARIMAUtils._
import common.ARIMAIdTransUtils._
import common.Innovations._
import stat.StatDenseVector
import stat.GradientApproximation._

class ARIMA(data: Array[Double], order: (Int, Int, Int), method: String, xv: List[String] = List(), pAppro: Int = 30, demean: Boolean = true,
  MeanMLEQ: Boolean = false, MaxLag: Int = 30, start_params: DenseVector[Double] = DenseVector.zeros[Double](1)) extends Preprocessing(data: Array[Double]) {

  private val p = order._1
  private val d = order._2
  private val q = order._3
  private val k = max(p, q)
  private val m = 20 + p + q
  private val initphi = AutoRegression.train(data, m).phi
  private val z = DenseVector.zeros[Double](n)

  def checkOrder: Boolean = p + q + 10 > n - 1 || p < 0 || q < 0 || p > maxOrder || q > maxOrder

  def Hannan_Rissanen(y : DenseVector[Double]): (DenseVector[Double], DenseVector[Double], Double, DenseVector[Double], DenseVector[Double]) = {
    val n = y.length
    val Z = DenseMatrix.zeros[Double](n - m - k, p + q)
    def slice(i: Int, x: DenseVector[Double], pq: Int): DenseVector[Double] = reverse(x((m + k + i - pq) to (m + k + i - 1)))
    
    for (i <- m to n - 1) z(i) = ARErrEst(y, i, initphi)
    for (i <- 0 to n - m - k - 1) {
      val temp = slice(i, z, q)
      for (j <- p to p + q - 1) Z(i, j) = temp(j - p)
    }
    if (p > 0) {
      for (i <- 0 to n - m - k - 1) {
        val temp = slice(i, y, p)
        for (j <- 0 to p - 1) Z(i, j) = temp(j)
      }
    }
    val G = inv(Z.t * Z)
    val b = G * Z.t * y((m + k) to n - 1)
    val xhat = Z * b
    val err = y((m + k) to n - 1) - xhat
    val sigma2 = sum(err.map(e => pow(e, 2))) / (n - m - k)
    val se = sqrt(diag(G) * sigma2)
    val phi = b(0 to p - 1)
    val se_phi = se(0 to p - 1)
    val theta = b(p to p + q - 1)
    val se_theta = se(p to p + q - 1)
    (phi, theta, sigma2, se_phi, se_theta)
  }

  def FitARMA(z: DenseVector[Double], p: Int, q: Int, pAppro: Int = 30, demean: Boolean = true,
    MeanMLEQ: Boolean = false, MaxLag: Int = 30, start_params: DenseVector[Double] = DenseVector.zeros[Double](1))
  : (Double, DenseVector[Double], DenseVector[Double], Double, Double, Double, DenseVector[Double], DenseMatrix[Double], DenseVector[Double], DenseVector[Double], Boolean, Int, Boolean, Boolean, (Int, Int, Int)) = {
    var Z = z
    var mz = mean(Z)
//    if (d > 0) Z = diff(x = z, lag = diffLag, differences = d)
    if (!demean) mz = 0
    val y = Z - mz
    var pApp = pAppro
    if (q == 0) pApp = p
    var ans = GetFitARMA(y = y, p = p, q = q, pAppro = pApp)
    var LL = ans._1
    var iter = 0
    var mu = 0.0
    var phiHat = ans._2
    var thetaHat = ans._3
    var convergence = ans._5
    if (MeanMLEQ && (p > 0 || q > 0)) {
      val MaxIter = 10
      var etol = 10.0
      while (etol > 1e-06 && iter < MaxIter) {
        val LLPrev = LL
        iter += 1
        phiHat = ans._2
        thetaHat = ans._3
        val g = TacvfARMA(phi = phiHat, theta = thetaHat, maxlag = pApp)
        val coefAR = PacfDL(c = g, LinearPredictor = true).right.get._2
        mu = GetARMeanMLE(y, coefAR)
        ans = GetFitARMA(y = y - mu, p = p, q = q, pAppro = pApp)
        LL = ans._1
        etol = abs(LL - LLPrev) / LLPrev
        convergence = ans._5
        println(convergence)
        if (!convergence) sys.error("GetARFit returned convergence = " + convergence)
      }
    }
    val muHat = mu + mz
    phiHat = ans._2
    thetaHat = ans._3
    var res = Z - muHat
    if (p > 0 || q > 0) {
      val g = TacvfARMA(phi = phiHat, theta = thetaHat, maxlag = pApp)
      val coefAR = PacfDL(c = g, LinearPredictor = true).right.get._2
      res = BackcastResidualsAR(y = y, phi = coefAR, Q = 100, demean = false)
    }
    val fits = Z - res
    val n = res.length
    val sigsq = sum(pow(res, 2)) / n
    var covHat = DenseVector.zeros[Double](0)
    // covariance matrix of the coefficient estimates
    // to be implemented
    // covHat = (InformationMatrixARMA(phiHat, thetaHat))/n}
    val resacf = new StatDenseVector(res)
    val racf = (resacf.acf(maxlag = MaxLag))(1 to MaxLag)
    val LBQ = LjungBoxTest(res = res, k = p + q, maxlag = MaxLag)
    val loglikelihood = ans._1
    (loglikelihood, phiHat, -thetaHat, sigsq, NaN, muHat, racf, LBQ, res, fits, demean, iter, convergence, MeanMLEQ, order)
  }

  def GetFitARMA(y: DenseVector[Double], p: Int, q: Int, pAppro: Int = 30,
    start_params: DenseVector[Double] = DenseVector.zeros[Double](1)): (Double, DenseVector[Double], DenseVector[Double], String, Boolean) = {
    var xinit = DenseVector.zeros[Double](p + q)
    var pApp = pAppro
    var phiHat = DenseVector.zeros[Double](0)
    var thetaHat = DenseVector.zeros[Double](0)
    var convergence = false
    var loglikelihood = 10000.0
    if (all(start_params)) {
      xinit = start_params
    }
    if (start_params.length != p + q) {
      println("GetARMAFit: init length not correct. Set to zero.")
      xinit = DenseVector.zeros[Double](p + q)
    }
    if (max(abs(start_params)) > 0.99) {
      println("GetARMAFit: init parameter setting outside (-0.99,0.99). Reset to 0.")
      xinit = DenseVector.zeros[Double](p + q)
    }
    val n = y.length
    val penaltyLoglikelihood = (-n / 2 * log(sum(pow(y, 2)) / n)) - 10000
    if (q == 0) {
      pApp = p
    }
    val CD = ChampernowneD(z = y, p = pApp, MeanZero = true)
    val LBFGS = new LBFGS[DenseVector[Double]](maxIter = 100, m = 7, tolerance = 1.0E-10)
    if (p > 0 && q > 0) {
      def EntropyARMA(x: DenseVector[Double]): Double = {
        if (max(abs(x)) > 0.99) -penaltyLoglikelihood + pow(max(abs(x)), 2)
        else {
          val zetaPhi = x(0 to p - 1)
          val zetaTheta = x(p to p + q - 1)
          val g = TacvfARMA(phi = PacfToAR(zeta = zetaPhi), theta = PacfToAR(zeta = zetaTheta), maxlag = pApp)
          val xpar = PacfDL(c = g, LinearPredictor = true).right.get._2
          -FastLoglikelihoodAR(phi = xpar, n = n, CD = CD)
        }
      }
      val gradient = new ApproximateGradientFunction(f = EntropyARMA, epsilon = 1.0E-10)
/*      val entropyFun = new DiffFunction[DenseVector[Double]] {
        def calculate(x: DenseVector[Double]) = {
          println(x)
          (EntropyARMA(x), gradientCenter4Approximate(EntropyARMA, x, 1.0E-5))
        }
      }*/
      var zpar = xinit
      var entropy = EntropyARMA(zpar)
      val entropyFundefault = new DiffFunction[DenseVector[Double]] {     
        def calculate(x: DenseVector[Double]) = {
           if (entropy > gradient.valueAt(x)) {
             zpar = x
             entropy = gradient.valueAt(x)
           }
          gradient.calculate(x)
        }
      }

      val statedef = LBFGS.minimizeAndReturnState(f = entropyFundefault, init = xinit)    
      // another optimizer could be used here, to be implmented.
      
      loglikelihood = entropy
      convergence = statedef.converged
      phiHat = PacfToAR(zpar(0 until p))
      thetaHat = PacfToAR(zpar(p until p + q))
    }
    
    if (p > 0 && q == 0) {
      def EntropyARMA(x: DenseVector[Double]): Double = {
        if (max(abs(x)) > 0.99)
          -penaltyLoglikelihood
        else {
          val zetaPhi = x(0 until p)
          val xpar = PacfToAR(zetaPhi)
          -FastLoglikelihoodAR(xpar, n, CD)
        }
      }
      val gradient = new ApproximateGradientFunction(f = EntropyARMA, epsilon = 1.0E-10)
      var zpar = xinit
      var entropy = EntropyARMA(zpar)
      val entropyFundefault = new DiffFunction[DenseVector[Double]] {     
        def calculate(x: DenseVector[Double]) = {
           if (entropy > gradient.valueAt(x)) {
             zpar = x
             entropy = gradient.valueAt(x)
           }
          gradient.calculate(x)
        }
      }

      val statedef = LBFGS.minimizeAndReturnState(f = entropyFundefault, init = xinit)    
      // another optimizer could be used here, to be implmented.
      
      loglikelihood = entropy
      convergence = statedef.converged
      phiHat = PacfToAR(zpar(0 until p))
      thetaHat = PacfToAR(zpar(p until p + q))
    }
    if (p == 0 && q > 0) {
      def EntropyARMA(x: DenseVector[Double]): Double = {
        if (max(abs(x)) > 0.99)
          -penaltyLoglikelihood
        else {
          val zetaTheta = x(0 to q - 1)
          val g = TacvfARMA(phi = DenseVector.zeros[Double](0), theta = PacfToAR(zeta = zetaTheta), maxlag = pApp)
          val xpar = PacfDL(c = g, LinearPredictor = true).right.get._2
          -FastLoglikelihoodAR(xpar, n, CD)
        }
      }
      val gradient = new ApproximateGradientFunction(f = EntropyARMA, epsilon = 1.0E-10)
      var zpar = xinit
      var entropy = EntropyARMA(zpar)
      val entropyFundefault = new DiffFunction[DenseVector[Double]] {     
        def calculate(x: DenseVector[Double]) = {
           if (entropy > gradient.valueAt(x)) {
             zpar = x
             entropy = gradient.valueAt(x)
           }
          gradient.calculate(x)
        }
      }

      val statedef = LBFGS.minimizeAndReturnState(f = entropyFundefault, init = xinit)    
      // another optimizer could be used here, to be implmented.
      
      loglikelihood = entropy
      convergence = statedef.converged
      phiHat = PacfToAR(zpar(0 until p))
      thetaHat = PacfToAR(zpar(p until p + q))
      //      println(phiHat)
      //      println(thetaHat)
    }
    (-loglikelihood, phiHat, thetaHat, "L-BFGS", convergence)
  }

  def fit_start_params(y : DenseVector[Double]): DenseVector[Double] = {
    var start_params_New = DenseVector.zeros[Double](p + q)
    if (start_params.length == 1 && start_params(0) == 0.0) {
      if (q != 0) {
        if (p != 0) {
          val (phi, theta, sigma2, se_phi, se_theta) = Hannan_Rissanen(y)
          start_params_New = DenseVector.vertcat(phi, theta)
        } else {
          val inittheta = MovingAverage.train(data = data, q = q, recursion_level = q).theta
          start_params_New = inittheta
        }
      } else {
        val initphi2 = AutoRegression.train(data = data, p = p).phi
        start_params_New = initphi2
      }
    } else {
      start_params_New = start_params
    }
    start_params_New
  }

  def run(): ARIMA_Model = {
    if (checkOrder) sys.error("Order p + q + 20 is too large, larger than max order n-1" +
      ",or q less than 1 , or larger than  min(n-1, 20*log10(n))")
    var resid = Resid(x, d, demean, xv)
    if (!demean) resid = Resid(x, d, demean, xv)
//    println(resid)
    val (phi, theta, sigma2, se_phi, se_theta) = Hannan_Rissanen(resid)
    val start_params = fit_start_params(resid)
    if (method.contains("hannan")) {
      val result_Hannan = (phi, theta, sigma2)
      val (new_sigma2, aicc) = innovation_update(x = resid, model = result_Hannan)
      val new_result_Hannan = (phi, theta, new_sigma2, aicc)
      new ARIMA_Model(x, xv, new_result_Hannan, order)
    } else {
      var result_MLE = FitARMA(resid, p, q, pAppro = pAppro, demean = demean,
        MeanMLEQ = MeanMLEQ, MaxLag = MaxLag, start_params = start_params)
      val model = (result_MLE._2, result_MLE._3,result_MLE._4)
      val (new_sigma2, aicc) = innovation_update(x = resid, model = model)
      result_MLE = result_MLE.copy(_4 = new_sigma2, _5 = aicc)
      new ARIMA_Model(x, xv, result_MLE)
    }
  }

}


class ARIMA_Model(y: DenseVector[Double], xv: List[String] = List(), 
    result: (Double, DenseVector[Double], DenseVector[Double], Double, Double, Double, DenseVector[Double], DenseMatrix[Double], DenseVector[Double], DenseVector[Double], Boolean, Int, Boolean, Boolean, (Int, Int, Int))) {

  def this(y: DenseVector[Double], xv: List[String], result: (DenseVector[Double], DenseVector[Double], Double, Double), order: (Int, Int, Int)) {
    this(y, xv, (0.0, result._1, result._2, result._3, result._4 , 0.0, DenseVector.zeros[Double](0), DenseMatrix.zeros[Double](0, 0),
      DenseVector.zeros[Double](0), DenseVector.zeros[Double](0), true, 0, false, false, order))
  }

  val (loglikelihood, phiHat, thetaHat, sigsq, aicc, muHat, racf, lbq, res, fits, demean, iter, convergence, meanMLEQ, order) =
    (result._1, result._2, result._3, result._4, result._5, result._6, result._7, result._8, result._9, result._10, result._11, result._12, result._13, result._14, result._15)

  val model = (phiHat, thetaHat, sigsq)
  val x = y
  val d = order._2

  def adf_Test() = {
    println("ADF Test = " + adfTest(y))
  }
  
  private def forecast_arma(y: DenseVector[Double], model: (DenseVector[Double], DenseVector[Double], Double), h: Int): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
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

  private def forecast_diff(x: DenseVector[Double], d: Int, xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = x.length
    var lag = 1
    if (xv(k - 1) == "diff") lag = xv(k).toInt
    val differences = d
    val y = diff(x, lag, differences)
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, h = h, k = k + 2)
    pred = diffinv(pred, lag, differences, x(n - lag*differences to n - 1))
    pred = pred(lag to lag + h - 1)
    if (phi == null) sys.error("diff before log")
    phi = DenseVector[Double](1.0 +: (-phi).toArray)
    phi = DenseVector.vertcat(phi, DenseVector.fill[Double](lag, 0.0)) - DenseVector.vertcat(DenseVector.fill[Double](lag, 0.0), phi)
    phi = -phi(1 until phi.length)
    (pred, phi, l, u)
  }

  
  private def forecast_log(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    var (pred, phi, l, u) = forecast_transform(x = log(x), xv = xv, model = model, h = h, k = k + 1)
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

  private def forecast_season(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int = 1): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
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
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, h = h, k = k + 2)
    pred = pred + s(n to n + h - 1)
    if (phi == null) {
      l = l + s(n to n + h - 1)
      u = u + s(n to n + h - 1)
    }
    (pred, phi, l, u)
  }

  private def forecast_trend(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int = 1): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val n = x.length
    val p = xv(k).toInt
    val X = DenseMatrix.zeros[Double](n + h, p + 1)
    for (j <- 0 to p) X(::, j) := pow(DenseVector.rangeD(1.0, n + h + 1), j)
    val b = X(0 until n, ::) \ x
    val xhat = X * b
    val y = x - xhat(0 until n)
    var (pred, phi, l, u) = forecast_transform(x = y, xv = xv, model = model, h = h, k = k + 2)
    pred = pred + xhat(n until n + h)
    if (phi == null) {
      l = l + xhat(n until n + h)
      u = u + xhat(n until n + h)
    }
    (pred, phi, l, u)
  }

  private def forecast_transform(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int = 1): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    if (k > xv.length) return forecast_arma(y = x, model = model, h = h)
    if (xv(k - 1) == "diff") return forecast_diff(x = x, d = d, xv = xv, model = model, h = h, k = k)
    if (xv(k - 1) == "log") return forecast_log(x = x, xv = xv, model = model, h = h, k = k)
    if (xv(k - 1) == "season") return forecast_season(x = x, xv = xv, model = model, h = h, k = k)
    if (xv(k - 1) == "trend") return forecast_trend(x = x, xv = xv, model = model, h = h, k = k)
    else sys.error("x vector transformation is invalid")
  }

  private def forecast(x: DenseVector[Double], xv: List[String] = List(), model: (DenseVector[Double], DenseVector[Double], Double), h: Int, k: Int = 1): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    var (pred, phi, l, u) = forecast_transform(x = x, xv = xv, model = model, h = h, k = 1)
    val se = DenseVector.zeros[Double](h)
    if (phi != null) {
      val theta = model._2
      val psi = ma_inf(phi = phi, theta = theta, n = h)
      val sigma2 = model._3
      def g(j: Int): Double = sum(pow(psi(0 until j), 2))
      for (i <- 1 to h) se(i - 1) = sqrt(sigma2 * g(i))
      l = pred - (se * 1.96)
      u = pred + (se * 1.96)
    }
    else return (pred, null, l, u)
    (pred, se, l, u)
  }

  def predict(predictionLength: Int = 12): (Array[Double], Array[Double], Array[Double], Array[Double]) = {
    val (pred, se, l, u) = forecast(x = y, xv = xv, model = model, h = predictionLength, k = 1)
    (pred.toArray, se.toArray, l.toArray, u.toArray)
  }
  
  def predict_DenseVector(predictionLength: Int = 12): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    val (pred, se, l, u) = forecast(x = y, xv = xv, model = model, h = predictionLength, k = 1)
    (pred, se, l, u)
  }
  
  def print_Predictions(predictionLength: Int = 12): Unit = {
    val (pred, se, l, u) = forecast(x = y, xv = xv, model = model, h = predictionLength, k = 1)
    println("length of predictions =" + predictionLength)
    println("predictions : " + pred)
    if (se != null) println("sqrt(MSE) = " + se)
    println("Lower Bound = " + l)
    println("Upper Bound = " + u)
  }

  def print_ARIMA_Model(): Unit = {
    val LL = loglikelihood
    var k = order._1 + order._3
    if (demean) k += 1
    val n = res.length
    val aic = -2 * LL + 2 * k
    val bic = -2 * LL + log(n) * k
    println("length of series = " + n + ",  number of parameters = " + k)
    println("differences = " + order._2)
    println("phi = " + phiHat)
    println("theta = " + thetaHat)
    println("sigma^2 = " + sigsq + ", mu = " + muHat)
    println("loglikelihood = " + LL + ",  aic = " + aic + ",  bic = " + bic + ", aicc = " + aicc)
  }
  
  def print_LjungBox(): Unit = {
    println("Ljung-Box : " + lbq)
  }
}

object ARIMA {
  def train(data: Array[Double], order: (Int, Int, Int) = (1, 0, 1), method: String = "MLE", xv: List[String] = List(), pAppro: Int = 30, demean: Boolean = true,
      MeanMLEQ: Boolean = false, MaxLag: Int = 30, start_params: DenseVector[Double] = DenseVector.zeros[Double](1)): ARIMA_Model = {
    new ARIMA(data = data, order = order, method = method, xv = xv, pAppro = pAppro, demean = demean, MeanMLEQ = MeanMLEQ, MaxLag = MaxLag, start_params).run()
  }
}