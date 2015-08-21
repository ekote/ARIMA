/**
 * @author  Jian Wang
 * @version 1.0
 * @date    10/06/2015
 */

package ARIMA

import IO.ParseJson._
import java.io._
import IO.DateFrequency._

// Test object for parsing Json and do forecasting
// the method is used to test and compare the results in my report
object TestJson {

  def main(args: Array[String]) {
    val data = parseJson(path = "D://sql/NCEORY.json", 
        startdate = "2011-01-10", enddate = "2013-09-10")
    val x = data.map(v => v._2)
    val x_weekly = FromDaily(x, 7)
    val x_monthly = FromDaily(x, 30)
    //    val model = ARIMA.train(x, order = (11,1,7), method="MLE", xv=List("diff","7","season","365"), demean =true)
    val model = AutoFit.autofit(x_monthly, pMin = 0, qMin = 0, pMax = 5, qMax = 5,
        xv = List("season", "12"), dLag = 1, dMin = 1, dMax = 1, Box_Cox = true)
    model.print_ARIMA_Model
    
    val forecast_array = model.predict(predictionLength = 8, enableBound = true,
        leastBound = 0, upperBound = 100000)
    val forecast_DenseVector = model.predict_DenseVector(predictionLength = 8, 
        enableBound = true, leastBound = 0, upperBound = 100000)    
    
    model.print_Predictions(8, enableBound = true)
    model.print_res_and_fits
    println(model.lbq)
    model.print_adf_Test
    model.print_train_data
    
    
    //    val model = AutoFit.autofit(x, pMin = 0, qMin = 0, pMax = 11, qMax = 11, xv = List("season", "365"), dLag = 7, dMin = 1, dMax = 1, Box_Cox = true)
//    val model = AutoFit.autofit(x, pMin = 0, qMin = 0, pMax = 10, qMax = 10, 
//        xv = List("season", "52"), dLag = 1, dMin = 1, dMax = 1, Box_Cox = true)
  }
  
}