
package ARIMA

object AutoFit {

  def autofit(x: Array[Double], xv: List[String] = List(), dMax: Int = 2, pMin: Int = 0, pMax: Int = 5, qMin: Int = 0, qMax: Int = 5): ARIMA_Model = {
    var model = ARIMA.train(x, order = (1, 0, 0), xv = List(), method = "MLE", demean = true)
    var maxAICc = model.aicc
    for (d <- 0 to dMax) {
      for (p <- pMin to pMax) {
        for (q <- qMin to qMax) {
          if (p != 0 || q != 0) {
            if (d == 0) {
              val model2 = ARIMA.train(x, order = (p, 0, q), xv = xv, method = "MLE", demean = true)
              if (model2.aicc < maxAICc) {
                model = model2
                maxAICc = model.aicc
              }
            } else {
              val model2 = ARIMA.train(x, order = (p, d, q), xv = List("diff", "1") ::: xv, method = "MLE", demean = true)
              if (model2.aicc < maxAICc){
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