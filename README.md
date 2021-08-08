# FVM
## Método de Volúmenes Finitos 

Código de Simulación del modelo bidominio y monodominio, respectivamente, a partir del Método de Volúmenes Finitos (FVM) con desarrollo en el entorno Matlab-R2019.

El FVM e integración explícita de Euler se emplean en el presente proyecto para discretizar la dinámica espacio-temporal y resolver numéricamente tanto las ecuaciones del bidominio como del monodominio, respectivamente. Al final se presenta los gráficos resultantes de simulaciones numéricas realizadas con código desarrollado en el ambiente de trabajo Matlab® para ambos modelos y diseñado en el esquema de volúmenes finitos.

Una importante propiedad de FVM es que los principios de conservación de masa, momento y energía, son preservados por las ecuaciones discretas deducidas de su aplicación, es decir, para el tema de estudio de la propagación de la actividad eléctrica cardíaca, el flujo se conserva con el método de volúmenes finitos. De forma general el FVM cimienta su base en el siguiente algoritmo de trabajo:

  - Descomponer el dominio en volúmenes de control, la malla admisible.
  - Formular las ecuaciones integrales de conservación para cada volumen de control.
  - Aproximar numéricamente las integrales.
  - Aproximar los valores de las variables en las caras y derivadas con la información de las variables nodales.
  - Ensamblar y resolver el sistema algebraico obtenido.

## Captura de Pantalla - Simulación - Método de Volúmenes Finitos

<div align="center">
    <center>
        <img src="simulate/figura35.jpg" width="700">
    </center>
</div>
<br>

