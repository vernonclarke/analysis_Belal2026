<h1 align="center">Curve Fitting Equations</h1>

## Definitions and Formulae

## <center>These routines are designed to provide fits for product functions</center>


### product function:

$$
\boldsymbol{ y = A (1 - e^{-t/\tau_1}) e^{-t/\tau_2} }
$$

This equation can be written in the alternative form:

$$
\boldsymbol{ y = A  (e^{-t/\tau_{decay}} - e^{-t/\tau_{rise}}) }
$$


where 

$$
\boldsymbol{ \tau_{rise} = \tau_1 \tau_2 / (\tau_1 + \tau_2) }
$$

$$
\boldsymbol{ \tau_{decay} = \tau_2 }
$$


### sum of two product functions:


$$
y = A_1 (1 - e^{-t/\tau_1}) e^{-t/\tau_2} + A_2 (1 - e^{-t/\tau_3}) e^{-t/\tau_4} 
$$

This equation can be written in the alternative form:

$$
y = A_1  (e^{-t/\tau_{decay1}} - e^{-t/\tau_{rise1}}) + A_2  (e^{-t/\tau_{decay2}} - e^{-t/\tau_{rise2}}) 
$$

where     

$$
\tau_{rise1} = \frac{\tau_1 \tau_2}{\tau_1 + \tau_2} 
$$

$$
\tau_{rise2} = \frac{\tau_3 \tau_4}{\tau_3 + \tau_4}
$$

$$
\tau_{decay1} = \tau_2 
$$

$$
\tau_{decay2} = \tau_4
$$


### alpha function:

$$
y = A t e^{-t/\tau} 
$$

### sum of two alpha functions:

$$
y = A_1 t e^{-t/\tau_1} + A_2 t e^{-t/\tau_2} 
$$

## <center>Solutions for the product function</center>

### Product function takes the form:

$$
y = A (1 - e^{-t/\tau_1}) e^{-t/\tau_2} 
$$

This equation can be written:

$$
y = A  (e^{-t/\tau_{decay}} - e^{-t/\tau_{rise}}) 
$$

where    
 
$$
\tau_{rise} = \tau_1 \tau_2 / (\tau_1 + \tau_2)
$$

$$
\tau_{decay} = \tau_2
$$

### Time to peak of response for product function:

In order to calculate $t_{peak}$, differentiate y with respect to t:

$$
\frac{dy}{dt} = A \left( \frac{e^{-t/\tau_{rise}}}{\tau_{rise}} - \frac{e^{-t/\tau_{decay}}}{\tau_{decay}} \right) 
$$

The time of the peak of the response $t = t_{peak}$ can be found by solving $\frac{dy}{dt} = 0$:

$$
0 = A \left( \frac{e^{-t_{peak}/\tau_{rise}}}{\tau_{rise}} - \frac{e^{-t_{peak}/\tau_{decay}}}{\tau_{decay}} \right)
$$

simplifying:

$$
\frac{e^{-t_{peak}/\tau_{decay}}}{\tau_{decay}} = \frac{e^{-t_{peak}/\tau_{rise}}}{\tau_{rise}}
$$

cross-multiplying:

$$
\frac{\tau_{decay}}{\tau_{rise}} = \frac{e^{-t_{peak}/\tau_{decay}}}{e^{-t_{peak}/\tau_{rise}}}
$$

simplifying:

$$
\frac{\tau_{decay}}{\tau_{rise}} = {e^{t_{peak} \left(\frac{\tau_{decay} - \tau_{rise}}{\tau_{decay} \cdot \tau_{rise}}\right)}}
$$

taking the natural logarithm of both sides and rearranging gives an expression for the time to peak $t_{peak}$:

$$
t_{peak} = \frac{\tau_{decay} \tau_{rise}}{\tau_{decay} - \tau_{rise}} ln\left(\frac{\tau_{decay}}{\tau_{rise}}\right)
$$

substituting for $\tau_{rise}$ and $\tau_{decay}$ gives an equivalent form in terms of $\tau_1$ and $\tau_2$:

$$
t_{peak} = \tau_1 ln\left(\frac{\tau_1 + \tau_2}{\tau_1}\right)
$$

To find the peak of the response $A_{peak}$, find solution where $t = t_{peak}$

$$
\boldsymbol{ A_{peak} = Af }
$$

where the fraction f is given by 

$$
f = {e^{-t_{peak}/\tau_{decay}} - e^{-t_{peak}/\tau_{rise}}}
$$

f in terms of $\tau_{decay}$ and $\tau_{rise}$: 

$$
f = {e^{-\frac{\tau_{rise}}{\tau_{decay} - \tau_{rise}} ln\left(\frac{\tau_{decay}}{\tau_{rise}}\right)} - e^{-\frac{\tau_{decay} }{\tau_{decay} - \tau_{rise}} ln\left(\frac{\tau_{decay}}{\tau_{rise}}\right)}}
$$

since $e^{-xln(y)} = e^{ln(y^{-x})} = y^{-x}$

$$
\boldsymbol{ f = {\left( \left(\frac{\tau_{decay}}{\tau_{rise}}\right)^{-\frac{\tau_{rise}}{\tau_{decay}-\tau_{rise}}} \right) \left( 1 - \frac{\tau_{rise}}{\tau_{decay}}\right) } }
$$

similarly in terms of $\tau_1$ and $\tau_2$: 

$$
\boldsymbol{ f = {\left( \left( \frac{\tau_1}{\tau_1+\tau_2} \right)^{\frac{\tau_1}{\tau_2}} \right) \frac{\tau_2}{\tau_1+\tau_2}} }
$$

### Area under the curve for the product function:

To find the area under the curve for the equation:

$$ 
y = A(e^{-t/\tau_{decay}} - e^{-t/\tau_{rise}}) 
$$

Integrate the function with respect to t then calculate the integral of y from 0 to $\infty$ (i.e. the area under the curve):

$$ 
\text{Area} = \int_{0}^{\infty} A(e^{-t/\tau_{decay}} - e^{-t/\tau_{rise}})dt 
$$

Solve this integral:

$$
\text{Area} = \int_{0}^{\infty} A(e^{-t/\tau_{decay}} - e^{-t/\tau_{rise}})dt = A \left[ -\tau_{decay} e^{-t/\tau_{decay}} + \tau_{rise} e^{-t/\tau_{rise}} \right]_0^{\infty} 
$$   

$$
= A (\tau_{decay} - \tau_{rise}) 
$$    

Area under the curve is given by

$$ 
\boldsymbol{ \text{Area} = A (\tau_{decay} - \tau_{rise}) = \frac{A_{peak}}{f} (\tau_{decay} - \tau_{rise}) }
$$

Similarly, in terms of $\tau_1$ and $\tau_2$:

$$ 
\boldsymbol{ \text{Area} = A \left(\frac{\tau_2^2}{\tau_1 + \tau_2}\right) = \frac{A_{peak}}{f} \left(\frac{\tau_2^2}{\tau_1 + \tau_2}\right) }
$$

where $A_{peak}$ is the peak amplitude and f is as previously defined (see above)

### Rise and decay kinetics:
Let p be the relative amplitudes of the response at some time t such that $p = y / A_{peak}$
rearranging the equation:

$$
y = \frac{A_{peak}}{f}(e^{-t/\tau_{decay}} - e^{-t/\tau_{rise}}) 
$$

gives:

$$
e^{-t/\tau_{decay}} - e^{-t/\tau_{rise}} -fp = 0
$$

This equation can be solved for t using a numerical root-finding algorithm (e.g. using an iterative method with an initialised variables as the starting value for the iteration):

Solving for $p_1$ and $p_2$ gives $t_1$ and $t_2$

if $p_2 > p_1$ then rise time is given by:

$$
rise_{p_1 - p_2} =  t_2 - t_1
$$

if $p_1 > p_2$ then decay time is given by:

$$
decay_{p_1 - p_2} =  t_2 - t_1
$$

for instance if $p_1$ = 0.2 and $p_2$ = 0.8 then $t_2$ - $t_1$ gives the 20 - 80 % rise time

likewise if $p_1$ = 0.9 and $p_2$ = 0.1 then $t_2$ - $t_1$ gives the 90 - 10 % decay time

## <center>Solutions for the alpha function (reference only; not implemented)</center>

### alpha function takes the form:

$$
y = A t e^{-t/\tau} 
$$

The solutions for the peak response ($A_{peak}$), time to peak ($t_{peak}$) and area are easier to calculate:

### Time to peak of response for alpha function:

In order to calculate $t_{peak}$, differentiate y with respect to t:

$$
\frac{dy}{dt} = A \left(1 - \frac{t}{\tau} \right) e^{-t/\tau} 
$$

The time of the peak of the response $t = t_{peak}$ can be found by solving $\frac{dy}{dt} = 0$:

$$
0 = A \left(1 - \frac{t}{\tau} \right) e^{-t/\tau} 
$$

simplifying:

$$
\boldsymbol{ t_{peak} = \tau }
$$

To find the peak of the response $A_{peak}$, find solution where $t = t_{peak}$

$$
\boldsymbol{ A_{peak} = A \tau e^{-1} }
$$

Substituting for $A = \frac{A_{peak}} {\tau} e^1$ in original equation gives an often used form of the alpha function 

$$
\boldsymbol{ y = A_{peak} \frac{t}{\tau} e^{1-t/\tau} }
$$

### Area under the curve for the alpha function:

To find the area under the curve for the equation:

$$ 
y = A t e^{-t/\tau} 
$$

Integrate the function with respect to t then calculate the integral of y from 0 to $\infty$ (i.e. the area under the curve):

$$ 
\text{Area} = \int_{0}^{\infty} A t e^{-t/\tau} dt 
$$

Solve this integral:

$$
\lim_{T \to \infty} \int_{0}^{T} Ate^{-t/\tau} dt = \lim_{T \to \infty} \left[ -A\tau e^{-t/\tau}(t+\tau) \right]_0^T 
$$ 

$$ 
= \lim_{T \to \infty} \left[ -A\tau e^{-T/\tau}(T+\tau) + A\tau^2e^{0} \right] 
$$

$$
\boldsymbol{ \text{Area} = A\tau^2 = A_{peak} \tau e^1 }
$$
    
### Rise and decay kinetics:
Let p be the relative amplitudes of the response at some time t such that $p = y / A_{peak}$
rearranging the equation:

$$
\boldsymbol{ y = A_{peak} \frac{t}{\tau} e^{1-t/\tau} }
$$

gives:

$$
 te^{-t/\tau} - p\tau e^{-1} = 0
$$

This equation can be solved for t using a numerical root-finding algorithm (e.g. using an iterative method with initialised variables as the starting values for the iteration):

Solving for $p_1$ and $p_2$ gives $t_1$ and $t_2$

$$
rise_{p_1 - p_2} =  t_2 - t_1
$$

if $p_1 > p_2$ then decay time is given by:

$$
decay_{p_1 - p_2} =  t_2 - t_1
$$

for instance if $p_1 = 0.2$ and $p_2 = 0.8$ then $t_2$ - $t_1$ gives the 20 - 80 % rise time

likewise if $p_1 = 0.9$ and $p_2 = 0.1$ then $t_2$ - $t_1$ gives the 90 - 10 % decay time






$A_1$, $τ_{rise}$, $τ_{decay}$, $t_{peak}$, delay, $r_{10-90}$, $d_{90-10}$, $area_1$
