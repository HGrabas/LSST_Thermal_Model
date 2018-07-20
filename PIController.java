package org.lsst.ccs.subsystem.rafts;

/**
 * Class: PI_Controller.java; Proportional, Integral controller <p> Description:
 * This PI_Controller class implements a proportional, integral feedback loop
 * controller for use in the thermal control of the LSST camera electronics. The
 * functionality includes anti-windup.
 *
 * @author Innes
 */
class PIController {

    // Parameters
    private double gain = 15;         // loop gain
    private double timeConst;         // integration time constant
    private double smoothTime = 100.; // input smoothing time in seconds   
    private double maxOutput = 130.0; // maximum PID output
    private double minOutput = 0.0;   // minimum PID output
    private double maxInput = 0.0;    // maximum input - limit setpoint to this
    private double minInput = 0.0;    // minimum input - limit setpoint to this
    private double tolerance = 0.1;   // error that is considered on target
    private double setpoint = 140.0;  // the desired average of the input
    private double awgain = 4.;       // gain in the anti-windup loop
    private double baseOutput = 0.;   // base output

    // Fields
    private double aveInput;	      // average input value
    private double smoothInput;       // smoothed input value
    private double measTime;          // time of measurement   
    private double lastTime = 0.;     // previous measurement time
    private double errorIntegral = 0.0; //sum of errors used in the integ. calc
    private double setError = 0.0;  // diff. between actual & desired values 


   /**
    * Allocate a PI object with the given constants for P, I
    *
    * @param Init_gain the loop gain coefficient
    * @param Init_timeConstant the integration time constant (seconds)
    */
    public PIController(double Init_gain, double Init_timeConstant) {

        gain = Init_gain;
        timeConst = Init_timeConstant;
    }


   /**
    * Return the current PI result. This is constrained by the max and min
    * outs
    *
    * @return the latest calculated output
    */
    public double performPI(double[] input, double time) {
        aveInput = 0.;
        for (double mInput : input) {
            aveInput += mInput / input.length;
        }
        measTime = time;
        return calculate();
    }


   /**
    * Set the PI Controller gain parameters.
    * Set the proportional, integral, and differential coefficients.
    *
    * @param p Proportional coefficient (overall loop gain)
    * @param intTime Integral coefficient (integration time constant)
    */
    public void setPID(double p, double intTime) {
        gain = p;
        timeConst = intTime;
    }


   /**
    * Get the Proportional coefficient
    *
    * @return proportional coefficient
    */
    public double getP() {
        return gain;
    }


   /**
    * Get the Integral coefficient
    *
    * @return integral coefficient
    */
    public double getI() {
        return timeConst;
    }


   /**
    * Sets the smoothing time.
    *
    * @param time the smoothing time
    */
    public void setSmoothTime(double time) {
        smoothTime = time;
    }


   /**
    * Sets the anti-windup gain.
    *
    * @param gain the anti-windup gain
    */
    public void setAwGain(double gain) {
        awgain = gain;
    }


   /**
    * Sets the base output.
    *
    * @param value the base output value
    */
    public void setBaseOutput(double value) {
        baseOutput = value;
    }


   /**
    * Set the percentage error which is considered tolerable for use with
    * OnTarget. (Input of 15.0 = 15 percent)
    *
    * @param percent error which is tolerable
    */
    public void setTolerance(double dK) {
        tolerance = dK;
    }


    /**
    * Sets the maximum and minimum values expected from the input.
    *
    * @param minimumInput the minimum value expected from the input
    * @param maximumInput the maximum value expected from the output
    */
    public void setInputRange(double minimumInput, double maximumInput) {
        minInput = minimumInput;
        maxInput = maximumInput;
        setSetpoint(setpoint);
    }


   /**
    * Sets the minimum and maximum values to write.
    *
    * @param minimumOutput the minimum value to write to the output
    * @param maximumOutput the maximum value to write to the output
    */
    public void setOutputRange(double minimumOutput, double maximumOutput) {
        minOutput = minimumOutput;
        maxOutput = maximumOutput;
    }


   /**
    * Set the setpoint for the PIDController
    *
    * @param setpoint the desired setpoint
    */
    public void setSetpoint(double setsetpoint) {
        setpoint = (maxInput <= minInput) ? setsetpoint :
                   (setsetpoint > maxInput) ? maxInput :
                   (setsetpoint < minInput) ? minInput : setsetpoint;
    }


   /**
    * Returns the current setpoint of the PIDController
    *
    * @return the current setpoint
    */
    public double getSetpoint() {
        return setpoint;
    }


   /**
    * Returns the current difference of the input from the setpoint
    *
    * @return the current error
    */
    public double getError() {
        return setError;
    }


   /**
    * Returns the current error integral
    *
    * @return the current error integral
    */
    public double getIntegral() {
        return errorIntegral;
    }


   /**
    * Sets the error integral
    *
    * @param value the error integral value to set
    */
    public void setIntegral(double value) {
        errorIntegral = value;
    }


   /**
    * Returns the current smoothed input
    *
    * @return the current smoothed input
    */
    public double getSmooth() {
        return smoothInput;
    }


   /**
    * Return true if the error is within the percentage of the input range,
    * determined by setTolerance.  This assumes that the maximum and minimum
    * input were set using setInput.
    *
    * @return true if the error is less than the tolerance
    */
    public boolean onTarget() {
        return (Math.abs(setError) < tolerance / 100 * (maxInput - minInput));
    }


   /**
    * Reset the previous error, the integral, and disable the controller.
    */
    public void reset() {
        errorIntegral = 0;
        lastTime = 0.;
    }


   /**
    * Smooth the input, perform the PI algorithm. This should only be called by
    * the performPI method
    */
    private double calculate() {
        // controller calculations
        // Calculate the error signal

        // apply a single pole RC filter to the input value

        if (lastTime == 0.) {
            lastTime = measTime - 1.;
            smoothInput = aveInput;
        }                            // initialize the filter

        // perform the filter
        smoothInput = ((measTime - lastTime) * aveInput
            + (smoothTime - measTime + lastTime) * smoothInput) / smoothTime;

        // find the "error"
        setError = (setpoint - smoothInput);

        // apply the overall loop gain
        double smoothOutput = gain * setError;

        // add the proportional and integral terms
        double propInt = smoothOutput + errorIntegral;

        // Make sure the final result is within bounds, truncate otherwise
        double output = (propInt > maxOutput) ? maxOutput :
                        (propInt < minOutput) ? minOutput : propInt;

        // calculate the amount of output truncation
        double outputTrunc = awgain * (output - propInt);

        // integrate the error less amount proportinal to the truncation
        errorIntegral +=
            (smoothOutput + outputTrunc) * (measTime - lastTime) / timeConst;

        lastTime = measTime;

        // add base output
        return output + baseOutput;
    }

}
