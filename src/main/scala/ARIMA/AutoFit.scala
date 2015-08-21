
/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

package ARIMA

// AutoFit object
object AutoFit {
/* autofit method to automatically fit ARIMA model
 * you have to determine the seasonality and trend before use it, as well as period and trend type (linear?)
 * then arguments:
 * x: the time series
 * xv: The data transformation list:
 * List("diff", 1, log, "season", 12, "trend",1) means take one difference with lag = 1, take log transformation (normally we use Box-Cox transformation)
 * and detect seasonality as period = 12 and trend component = 1 (linear).. plese follow the order as indicated, and just include what you want to perform or want
 * time series to be transformed
 * dMin: minimum difference
 * (suggestion: do figure out the order and the lag of difference in advance because for different differences, aicc is not comparable
 * that is set dMin=dMax=some number)
 * dMax: maximum difference
 * dLag: number of lag in difference
 * pMin: min order p
 * pMax: max order p
 * qMin: min order q
 * qMax: max order q
 * method: "hannan" or "MLE" (specify hannan otherwise MLE)
 * pAppro: the order AR model approximation of infinite, 30 default, according to your time series length, if 
 * your time series is short, change it to smaller number.
 * demean: if TRUE, mean parameter included otherwise assumed zero
 * MeanMLEQ: exact MLE for mean, ignored unless demean=true
 * Box_Cox: Use Box-Cox Transformation or not
 * Lambda: the lambda of box-cox transformation if Box_Cox=true  
 * MaxLag: maximum number of lags for portmanteau test
 * 
 * return:
 * a auto fitted ARIMA model with lowest AICc value
 */
  def autofit(x: Array[Double], xv: List[String] = List(), dMin: Int = 0, dMax: Int = 1, dLag: Int = 1, pMin: Int = 0, pMax: Int = 5, 
      qMin: Int = 0, qMax: Int = 5, method: String = "MLE", MeanMLEQ: Boolean = false, demean: Boolean = true, Box_Cox: Boolean = false, Lambda: Double = 0.5
      , MaxLag: Int = 30): ARIMA_Model = {
    var model = ARIMA.train(x, order = (pMin, dMin, qMin), xv = xv, method = method, demean = demean, MeanMLEQ = MeanMLEQ, Box_Cox = Box_Cox, Lambda = Lambda, MaxLag = MaxLag)
    var maxAICc = model.aicc
    for (d <- dMin to dMax) {
      for (p <- pMin to pMax) {
        for (q <- qMin to qMax) {
          if (p != 0 || q != 0) {
            if (d == 0) {
              val model2 = ARIMA.train(x, order = (p, 0, q), xv = xv, method = method, demean = demean, MeanMLEQ = MeanMLEQ, Box_Cox = Box_Cox, Lambda = Lambda,  MaxLag = MaxLag)
              if (model2.aicc < maxAICc) {
                model = model2
                maxAICc = model.aicc
              }
            } else {
              val model2 = ARIMA.train(x, order = (p, d, q), xv = List("diff", dLag + "") ::: xv, method = method, demean = demean, MeanMLEQ = MeanMLEQ, Box_Cox = Box_Cox, Lambda = Lambda, MaxLag = MaxLag)
              if (model2.aicc < maxAICc) {
                model = model2
                maxAICc = model.aicc
              }
            }
          }
        }
      }
    }
    model
  }

}