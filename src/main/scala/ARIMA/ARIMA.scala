
/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

// package TimeSeriesForecast.ARIMA

package ARIMA

import breeze.linalg.{ DenseVector, DenseMatrix, max, inv, diag, sum, all, reverse }
import breeze.numerics.{ abs, log, NaN, pow, sqrt, exp }
import breeze.optimize.{ ApproximateGradientFunction, LBFGS, DiffFunction }
import breeze.stats.mean
import common.ARIMAUtils._
import common.ARIMAIdTransUtils._
import common.Innovations._
import common.Forecast._
import stat.StatDenseVector
import stat.GradientApproximation._

/*object: ARIMA model
 * this is the main class for ARIMA model, Fits an ARIMA(p,d,q) model using the algorithm given in McLeod and Zhang (2007) or exact MLE
 * The input arguments:
 * 
 * data: the time series
 * order: p,d,q
 * method: "hannan" or "MLE" (specify hannan otherwise MLE)
 * xv: The data transformation list 
 * pAppro: the order AR model approximation of infinite, 30 default, according to your time series length, if 
 * your time series is short, change it to smaller number.
 * demean: if TRUE, mean parameter included otherwise assumed zero
 * MeanMLEQ: exact MLE for mean, ignored unless demean=true
 * Box_Cox: Use Box-Cox Transformation or not
 * Lambda: the lambda of box-cox transformation if Box_Cox=true  
 * MaxLag: maximum number of lags for portmanteau test
 * start_params: start parameters of estimation of parameters using MLE can be provided here otherwise default using 
 * Hannan to calculate start parameters
*/
class ARIMA(data: Array[Double], order: (Int, Int, Int), method: String, xv: List[String] = List(), pAppro: Int = 30, demean: Boolean = true,
  MeanMLEQ: Boolean = false, Box_Cox: Boolean = false, Lambda: Double = 0.5, MaxLag: Int = 30, start_params: DenseVector[Double] = DenseVector.zeros[Double](1)) extends Preprocessing(data: Array[Double]) {

  private val p = order._1
  private val d = order._2
  private val q = order._3
  private val k = max(p, q)
  private val m = 20 + p + q
  private val initphi = AutoRegression.train(data, m).phi
  private val z = DenseVector.zeros[Double](n)

  // check orders provided are valid
  def checkOrder: Boolean = p + q + 10 > n - 1 || p < 0 || q < 0 || p > maxOrder || q > maxOrder

  //Hannan Rissanen algorithm, just need to provide time series data
  def Hannan_Rissanen(y: DenseVector[Double]): (DenseVector[Double], DenseVector[Double], Double, DenseVector[Double], DenseVector[Double]) = {
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

  /* Main method to estimate the ARMA/ARIMA model coefficients
   * Fit ARMA/ARIMA using fast MLE algorithm
   * Fits an ARIMA(p,d,q) model using the algorithm given in McLeod and Zhang (2007).
   * arguments:
   * z: tranformed time series
   * p: order of AR
   * q: order of MA
   * demean: if TRUE, mean parameter included otherwise assumed zero
   * MeanMLEQ: exact MLE for mean, ignored unless demean=true
   * pAppro: order of approximation to be used
   * MaxLag: maximum number of lags for portmanteau test
   * start_params: start parameters of estimation of parameters using MLE can be provided here otherwise default using 
   * Hannan to calculate start parameters
   * return: results of fitting model 
   * loglikelihood: loglikelihood estimated by MLE and LBFGS algorithm
   * phiHat: estimated AR parameters
   * thetaHat: estimated MA parameters
   * sigsq: innovation variance estimate
   * aicc: return NaN, will be calculated later
   * muHat: estimate of the mean 
   * racf: residual autocorrelations
   * LBQ: table of Ljung-Box portmanteau test statistics
   * res: innovation residuals, same length as z
   * fits: fitted values, same length as z
   * demean: TRUE if mean estimated otherwise assumed zero
   * iter: number of iterations in mean mle estimation
   * convergence: true if LBFGS convergence, otherwise LBFGS not convergence
   * MeanMLQ: TRUE if mle for mean algorithm used
   * order: model orders
  */
  def FitARMA(z: DenseVector[Double], p: Int, q: Int, pAppro: Int = 30, demean: Boolean = true,
    MeanMLEQ: Boolean = false, MaxLag: Int = 30, start_params: DenseVector[Double] = DenseVector.zeros[Double](1)): (Double, DenseVector[Double], DenseVector[Double], Double, Double, Double, DenseVector[Double], DenseMatrix[Double], DenseVector[Double], DenseVector[Double], Boolean, Int, Boolean, Boolean, (Int, Int, Int)) = {
    var Z = z
    var mz = mean(Z)
    //    if (d > 0) Z = diff(x = z, lag = diffLag, differences = d)
    if (!demean) mz = 0
    val y = Z - mz
    var pApp = pAppro
    if (q == 0) pApp = p
    var ans = GetFitARMA(y = y, p = p, q = q, pAppro = pApp, start_params = start_params)
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
        //        println(convergence)
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
    // can be used to calculate the coefficients standard deviation
    // covHat = (InformationMatrixARMA(phiHat, thetaHat))/n}
    val resacf = new StatDenseVector(res)
    val racf = (resacf.acf(maxlag = MaxLag))(1 to MaxLag)
    val LBQ = LjungBoxTest(res = res, k = p + q, maxlag = MaxLag)
    val loglikelihood = ans._1
    (loglikelihood, phiHat, -thetaHat, sigsq, NaN, muHat, racf, LBQ, res, fits, demean, iter, convergence, MeanMLEQ, order)
  }
  
  /*
   * Fit ARMA(p,q) model to mean zero time series.
   * The algorithm of McLeod and Zhang (2007) is used.
   * LBGFS as an optimization algorithm is used here to minimize an equation
   * Approximation Gradiance Function can be used here, if more accuracy needed, you can use my optimized implementation of 
   * Approximation of gradiance function.
   * y: time series data
   * p: order of AR model
   * q: order of MA model
   * pAppro: AR approximation
   * start_params: start parameters of estimation of parameters using MLE can be provided here otherwise default using 
   * Hannan to calculate start parameters
   */
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
    if (start_params.length != p + q || start_params.length == 0) {
      println("GetARMAFit: init length not correct. Set to zero.")
      xinit = DenseVector.zeros[Double](p + q)
    } else {
      if (max(abs(start_params)) > 0.99) {
        println("GetARMAFit: init parameter setting outside (-0.99,0.99). Reset to 0.")
        xinit = DenseVector.zeros[Double](p + q)
      }
    }
    val n = y.length
    val penaltyLoglikelihood = (-n / 2 * log(sum(pow(y, 2)) / n)) - 10000
    if (q == 0) {
      pApp = p
    }
    val CD = ChampernowneD(z = y, p = pApp, MeanZero = true)
    val LBFGS = new LBFGS[DenseVector[Double]](maxIter = 100000, m = 7, tolerance = 1.0E-10)
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

  // to fit the start parameters.
  // if start_params not provided, use default, Hannan as start parameters
  // if p = 0 or q = 0, then we will use AR model or MA model to get the start parameters.
  def fit_start_params(y: DenseVector[Double]): DenseVector[Double] = {
    var start_params_New = DenseVector.zeros[Double](p + q)
    if (start_params.length == 1 && start_params(0) == 0.0) {
      if (q != 0) {
        if (p != 0) {
          val Hannan_params = try {
            Some(Hannan_Rissanen(y))
          } catch {
            case e: Exception => None
          }
          val (phi, theta, sigma2, se_phi, se_theta) = Hannan_params match {
            case Some(a) => (a._1, a._2, a._3, a._4, a._5)
            case None => (DenseVector[Double](), DenseVector[Double](), NaN, DenseVector[Double](), DenseVector[Double]())
          }
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

  // train ARIMA model,
  // check order first then fit Hannan 
  // determine which method to use and train model, then return an ARIMA_Model class to use for forecasting
  // run innovation update to use innovation algorithm to find variance as well as AICc
  def run(): ARIMA_Model = {
    if (checkOrder) sys.error("Order p + q + 20 is too large, larger than max order n-1" +
      ",or q less than 1 , or larger than  min(n-1, 20*log10(n))")
    var xt = x
    if (Box_Cox) xt = BoxCox(x, Lambda)
    var resid = Resid(xt, d, demean, xv)
    val Hannan_params = try {
      if (p != 0 || q != 0) Some(Hannan_Rissanen(resid))
      else Some((DenseVector[Double](), DenseVector[Double](), NaN, DenseVector[Double](), DenseVector[Double]()))
    } catch {
      case e: Exception => None
    }
    val (phi, theta, sigma2, se_phi, se_theta) = Hannan_params match {
      case Some(a) => (a._1, a._2, a._3, a._4, a._5)
      case None => (DenseVector[Double](), DenseVector[Double](), NaN, DenseVector[Double](), DenseVector[Double]())
    }
    if (method.contains("hannan")) {
      val result_Hannan = (phi, theta, sigma2)
      val (new_sigma2, aicc) = innovation_update(x = resid, model = result_Hannan)
      val new_result_Hannan = (phi, theta, new_sigma2, aicc)
      new ARIMA_Model(x, resid, xv, Box_Cox, Lambda, new_result_Hannan, order)
    } else {
      var start_params_ = start_params
      if (start_params.length == 1 && !all(start_params)) {
        start_params_ = fit_start_params(resid)
      }
      var result_MLE = FitARMA(resid, p, q, pAppro = pAppro, demean = demean,
        MeanMLEQ = MeanMLEQ, MaxLag = MaxLag, start_params = start_params_)
      val model = (result_MLE._2, result_MLE._3, result_MLE._4)
      val (new_sigma2, aicc) = innovation_update(x = resid, model = model)
      result_MLE = result_MLE.copy(_4 = new_sigma2, _5 = aicc)
      new ARIMA_Model(x, resid, xv, Box_Cox, Lambda, result_MLE)
    }
  }

}

/* Main class called by ARIMA_Model class, take model results as argument.
 * refer to the arguments in ARIMA class, and results of FitARMA methods (the return value)
 * you can perform forecast, predict, print model, ADF test and many other operations here
 * you can output results of model here
 */
class ARIMA_Model(y: DenseVector[Double], resid: DenseVector[Double], xv: List[String] = List(), Box_Cox: Boolean, Lambda: Double,
  result: (Double, DenseVector[Double], DenseVector[Double], Double, Double, Double, DenseVector[Double], DenseMatrix[Double], DenseVector[Double], DenseVector[Double], Boolean, Int, Boolean, Boolean, (Int, Int, Int))) {

  // Hannan algorithm initializaiton of ARIMA_Model class
  def this(y: DenseVector[Double], resid: DenseVector[Double], xv: List[String], Box_Cox: Boolean = false, Lambda: Double = 0.5, result: (DenseVector[Double], DenseVector[Double], Double, Double), order: (Int, Int, Int)) {
    this(y, resid, xv, Box_Cox, Lambda, (0.0, result._1, result._2, result._3, result._4, 0.0, DenseVector.zeros[Double](0), DenseMatrix.zeros[Double](0, 0),
      DenseVector.zeros[Double](0), DenseVector.zeros[Double](0), true, 0, false, false, order))
  }

  /* return the results of model
   * loglikelihood: loglikelihood estimated by MLE and LBFGS algorithm
   * phiHat: estimated AR parameters
   * thetaHat: estimated MA parameters
   * sigsq: innovation variance estimate
   * aicc: return aicc value, obtained by innovation update
   * muHat: estimate of the mean 
   * racf: residual autocorrelations
   * LBQ: table of Ljung-Box portmanteau test statistics
   * res: innovation residuals, same length as z
   * fits: fitted values, same length as z
   * demean: TRUE if mean estimated otherwise assumed zero
   * iter: number of iterations in mean mle estimation
   * convergence: true if LBFGS convergence, otherwise LBFGS not convergence
   * MeanMLQ: TRUE if mle for mean algorithm used
   * order: model orders
   * All of these value can be accessed globally, which means you can return them
   */
  val (loglikelihood, phiHat, thetaHat, sigsq, aicc, muHat, racf, lbq, res, fits, demean, iter, convergence, meanMLEQ, order) =
    (result._1, result._2, result._3, result._4, result._5, result._6, result._7, result._8, result._9, result._10, result._11, result._12, result._13, result._14, result._15)

  val model = (phiHat, thetaHat, sigsq)
  var xt = y
  if (Box_Cox) xt = BoxCox(y, Lambda)
  val d = order._2
  var xv2 = xv
  if (d != 0 && !xv.contains("diff")) xv2 = "diff" :: "1" :: xv

  // print adf test results
  def print_adf_Test() = {
    println("ADF Test = " + adfTest(resid))
  }
  // return adf test results
  // 
  def adf_Test() : (Double, Double, String, Double) = {
    adfTest(resid)
  }
  
  // print train data
  def print_train_data(): Unit = {
    println("Resids for constructing ARIMA model = " + resid)
  }
  
  /* predict with ARIMA model
   * predictionLenght: the length of prediction
   * enableBound: true if you want to set the boundary
   * leastBound: the lower bound of prediction
   * upperBound: the upperBound of prediction
   * return: Array[Double] 
   * prediction: predictions values
   * se: standard error
   * Lower Bound: the lower bound with 95% confidence interval
   * upper Bound: upper bound with 95% confidence interval
   */
  def predict(predictionLength: Int = 12, enableBound: Boolean = false, leastBound: Int = 0, upperBound: Int = 60000): (Array[Double], Array[Double], Array[Double], Array[Double]) = {
    var (pred, se, l, u) = forecast(x = xt, xv = xv2, model = model, d = d, h = predictionLength, k = 1, demean = demean, enableBound = enableBound, leastBound = leastBound, upperBound = upperBound)
    if (Box_Cox) {
      pred = InvBoxCox(pred, Lambda)
      l = InvBoxCox(l, Lambda)
      u = InvBoxCox(u, Lambda)
      se = (u - l) / (2 * 1.96)
    }
    if (se != null) (pred.toArray, se.toArray, l.toArray, u.toArray)
    else (pred.toArray, null, l.toArray, u.toArray)
  }

  // same as predict, just return DenseVector[Double] instead of Array[Doubel]
  def predict_DenseVector(predictionLength: Int = 12, enableBound: Boolean = false, leastBound: Int = 0, upperBound: Int = 60000): (DenseVector[Double], DenseVector[Double], DenseVector[Double], DenseVector[Double]) = {
    var (pred, se, l, u) = forecast(x = xt, xv = xv2, model = model, d = d, h = predictionLength, k = 1, demean = demean, enableBound = enableBound, leastBound = leastBound, upperBound = upperBound)
    if (Box_Cox) {
      pred = InvBoxCox(pred, Lambda)
      l = InvBoxCox(l, Lambda)
      u = InvBoxCox(u, Lambda)
      se = (u - l) / (2 * 1.96)
    }
    if (se != null) (pred, se, l, u)
    else (pred, null, l, u)
  }

  // same as predict_DenseVector() instead of return, just print the results
  def print_Predictions(predictionLength: Int = 12, enableBound: Boolean = false, leastBound: Int = 0, upperBound: Int = 60000): Unit = {
    var (pred, se, l, u) = forecast(x = xt, xv = xv2, model = model, d = d, h = predictionLength, k = 1, demean = demean, enableBound = enableBound, leastBound = leastBound, upperBound = upperBound)
    if (Box_Cox) {
      pred = InvBoxCox(pred, Lambda)
      l = InvBoxCox(l, Lambda)
      u = InvBoxCox(u, Lambda)
      se = (u - l) / (2 * 1.96)
    }
    println("length of predictions = " + predictionLength)
    println("predictions : " + pred)
    if (se != null) println("sqrt(MSE) = " + se)
    println("Lower Bound = " + l)
    println("Upper Bound = " + u)
  }

  // print ARIMA model,
  // length of series, number of parameters = p + q	
  // orders
  // phi values, theta values
  // sigma square variance and mean
  // loglikelihood, aic, bic, aicc
  def print_ARIMA_Model(): Unit = {
    val LL = loglikelihood
    var k = order._1 + order._3
    val numofparams = order._1 + order._3
    if (demean) k += 1
    val n = res.length
    val aic = -2 * LL + 2 * k
    val bic = -2 * LL + log(n) * k
    println("length of series = " + n + ",  number of parameters = " + numofparams)
    println("phi = " + order._1 + ", differences = " + order._2 + ", theta = " + order._3)
    println("phi = " + phiHat)
    println("theta = " + thetaHat)
    println("sigma^2 = " + sigsq + ", mu = " + muHat)
    println("loglikelihood = " + LL + ",  aic = " + aic + ",  bic = " + bic + ", aicc = " + aicc)
  }

  // print coefficients status
  // print phi, thata values, iteration times and convergence or not in LBFGS
  def print_coef_states(): Unit = {
    println("phi = " + phiHat)
    println("theta = " + thetaHat)
    println("iter = " + iter)
    println("convergence = " + convergence)
  }

  // print LjungBox result
  def print_LjungBox(): Unit = {
    println("Ljung-Box : " + lbq)
  }

  // print residuals, fits, and residual autocorrelations
  def print_res_and_fits(): Unit = {
    println("Residuals : " + res)
    println("fits : " + fits)
    println("residual autocorrelations : " + racf)
  }

  // return residuals and fits
  def model_res_and_fits(): (DenseVector[Double], DenseVector[Double]) = {
    (res, fits)
  }

}

/* build ARIMA object, and main method to be called. 
 * arugments:
 * data: the time series
 * order: p,d,q
 * method: "hannan" or "MLE" (specify hannan otherwise MLE)
 * xv: The data transformation list:
 * List("diff", 1, log, "season", 12, "trend",1) means take one difference with lag = 1, take log transformation (normally we use Box-Cox transformation)
 * and detect seasonality as period = 12 and trend component = 1 (linear).. plese follow the order as indicated, and just include what you want to perform or want
 * time series to be transformed
 * pAppro: the order AR model approximation of infinite, 30 default, according to your time series length, if 
 * your time series is short, change it to smaller number.
 * demean: if TRUE, mean parameter included otherwise assumed zero
 * MeanMLEQ: exact MLE for mean, ignored unless demean=true
 * Box_Cox: Use Box-Cox Transformation or not
 * Lambda: the lambda of box-cox transformation if Box_Cox=true  
 * MaxLag: maximum number of lags for portmanteau test
 * start_params: start parameters of estimation of parameters using MLE can be provided here otherwise default using 
 * Hannan to calculate start parameters
 */

object ARIMA {
  def train(data: Array[Double], order: (Int, Int, Int) = (1, 0, 1), method: String = "MLE", xv: List[String] = List(), pAppro: Int = 30, demean: Boolean = true,
    MeanMLEQ: Boolean = false, Box_Cox: Boolean = false, Lambda: Double = 0.5, MaxLag: Int = 30, start_params: DenseVector[Double] = DenseVector.zeros[Double](1)): ARIMA_Model = {
    new ARIMA(data = data, order = order, method = method, xv = xv, pAppro = pAppro, demean = demean, MeanMLEQ = MeanMLEQ, Box_Cox = Box_Cox, Lambda = Lambda, MaxLag = MaxLag, start_params).run()
  }
}