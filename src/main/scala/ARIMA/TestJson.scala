package ARIMA

import IO.ParseJson._
import java.io._
import IO.DateFrequency._

object TestJson {

  def main(args: Array[String]) {
    val data = parseJson(path = "D://sql/NCECDG.json", startdate = "2011-01-10",
        enddate = "2014-01-10")
    val x = data.map(v => v._2)
    val x_weekly = FromDaily(x, 7)
    val x_monthly = FromDaily(x, 30)
    //    val model = ARIMA.train(x, order = (11,1,7), method="MLE", xv=List("diff","7","season","365"), demean =true)
    val model = AutoFit.autofit(x_monthly, pMin = 0, qMin = 0, pMax = 5, qMax = 5, xv = List("season", "12"), dLag = 1, dMin = 1, dMax = 1, Box_Cox = true)
//    val model = AutoFit.autofit(x, pMin = 0, qMin = 0, pMax = 11, qMax = 11, xv = List("season", "365"), dLag = 7, dMin = 1, dMax = 1, Box_Cox = true)
//    val model = AutoFit.autofit(x, pMin = 0, qMin = 0, pMax = 10, qMax = 10, 
//        xv = List("season", "52"), dLag = 1, dMin = 1, dMax = 1, Box_Cox = true)
    model.print_ARIMA_Model
    model.print_Predictions(8, enableBound = true)
    model.print_res_and_fits
    model.adf_Test
    model.print_train_data
  }
  
}