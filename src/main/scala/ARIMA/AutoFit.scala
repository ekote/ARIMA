
package ARIMA

object AutoFit {

  def autofit(x: Array[Double], xv: List[String] = List(), dMin: Int = 0, dMax: Int = 1, dLag: Int = 1, pMin: Int = 0, pMax: Int = 5, 
      qMin: Int = 0, qMax: Int = 5, method: String = "MLE", MeanMLEQ: Boolean = false, demean: Boolean = true, Box_Cox: Boolean = false, Lambda: Double = 0.5): ARIMA_Model = {
    var model = ARIMA.train(x, order = (pMin, dMin, qMin), xv = xv, method = method, demean = demean, MeanMLEQ = MeanMLEQ, Box_Cox = Box_Cox, Lambda = Lambda)
    var maxAICc = model.aicc
    for (d <- dMin to dMax) {
      for (p <- pMin to pMax) {
        for (q <- qMin to qMax) {
          if (p != 0 || q != 0) {
            if (d == 0) {
              val model2 = ARIMA.train(x, order = (p, 0, q), xv = xv, method = method, demean = demean, MeanMLEQ = MeanMLEQ, Box_Cox = Box_Cox, Lambda = Lambda)
              if (model2.aicc < maxAICc) {
                model = model2
                maxAICc = model.aicc
              }
            } else {
              val model2 = ARIMA.train(x, order = (p, d, q), xv = List("diff", dLag + "") ::: xv, method = method, demean = demean, MeanMLEQ = MeanMLEQ, Box_Cox = Box_Cox, Lambda = Lambda)
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