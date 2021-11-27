import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from decimal import Decimal as dcm

# Две непересекающиеся окружности
def test1():
  theta = np.linspace(0, 2 * np.pi, 50)

  # Первая окружность
  r1 = np.sqrt(1)
  x1 = r1 * np.cos(theta) - 2
  y1 = r1 * np.sin(theta)

  # Вторая окружность
  r2 = np.sqrt(1)
  x2 = r2 * np.cos(theta) + 2
  y2 = r2 * np.sin(theta)

  a = np.arange(-2, 2, 0.1, dtype = np.float64)
  b = np.arange(-2, 2, 0.1, dtype = np.float64)

  xg, yg = np.meshgrid(a, b)

  z = ((xg + 2)**2 + yg**2 - 1)**2 + ((xg - 2)**2 + yg**2 - 1)**2

  xval = np.linspace(0, 200, 20)

  return x1, y1, x2, y2, xg, yg, z, xval

# Две окружности, пересекающиеся в одной точке
def test2():
  theta = np.linspace(0, 2 * np.pi, 50)

  # Первая окружность
  r1 = np.sqrt(1)
  x1 = r1 * np.cos(theta)
  y1 = r1 * np.sin(theta)

  # Вторая окружность
  r2 = np.sqrt(1)
  x2 = r2 * np.cos(theta) + 2
  y2 = r2 * np.sin(theta)

  a = np.arange(-1, 3, 0.01, dtype = np.float64)
  b = np.arange(-1, 1.5, 0.01, dtype = np.float64)

  xg, yg = np.meshgrid(a, b)

  z = (xg**2 + yg**2 - 1)**2 + ((xg - 2)**2 + yg**2 - 1)**2

  xval = np.linspace(0, 10, 20)

  return x1, y1, x2, y2, xg, yg, z, xval

# Две окружности, пересекающиеся в двух точках
def test3():
  theta = np.linspace(0, 2 * np.pi, 50)

  # Первая окружность
  r1 = np.sqrt(1)
  x1 = r1 * np.cos(theta) + 1
  y1 = r1 * np.sin(theta)

  # Вторая окружность
  r2 = np.sqrt(1)
  x2 = r2 * np.cos(theta) + 2
  y2 = r2 * np.sin(theta)

  a = np.arange(0, 3, 0.1, dtype = np.float64)
  b = np.arange(-1, 2, 0.1, dtype = np.float64)

  xg, yg = np.meshgrid(a, b)

  z = ((xg - 1)**2 + yg**2 - 1)**2 + ((xg - 2)**2 + yg**2 - 1)**2

  xval = np.linspace(0, 10, 20)

  return x1, y1, x2, y2, xg, yg, z, xval

# Две окружности, пересекающиеся в одной точке.
# К ним добавлена прямая
def test4():
  theta = np.linspace(0, 2 * np.pi, 50)

  # Первая окружность
  r1 = np.sqrt(1)
  x1 = r1 * np.cos(theta)
  y1 = r1 * np.sin(theta)

  # Вторая окружность
  r2 = np.sqrt(1)
  x2 = r2 * np.cos(theta) + 2
  y2 = r2 * np.sin(theta)

  xl = np.linspace(-0.5, 3, 10)
  yl = xl - 1

  a = np.arange(-1, 3, 0.01, dtype = np.float64)
  b = np.arange(-1.5, 1.5, 0.01, dtype = np.float64)

  xg, yg = np.meshgrid(a, b)

  z = (xg**2 + yg**2 - 1)**2 + ((xg - 2)**2 + yg**2 - 1)**2 + (yg - xg + 1)**2

  xval = np.linspace(0, 10, 20)

  return x1, y1, x2, y2, xl, yl, xg, yg, z, xval  

# Три попарно пересекающиеся прямые
def test5():
  x1 = np.linspace(-2, 2, 10)
  y1 = x1 + 1

  y2 = -x1

  y3 = -0.5 * x1

  a = np.arange(-2, 2, 0.01, dtype = np.float64)
  b = np.arange(-2, 3, 0.01, dtype = np.float64)

  xg, yg = np.meshgrid(a, b)

  z = (xg - yg + 1)**2 + (xg + yg)**2 + (0.5 * xg + yg)**2

  xval = np.linspace(0, 10, 20)

  return x1, y1, y2, y3, xg, yg, z, xval  

# Синусоида и прямая
def test6():
  x1 = np.linspace(-5, 5, 100)
  y1 = np.sin(x1)

  x2 = np.linspace(-3, 3, 10)
  y2 = x2

  a = np.arange(-4, 4, 0.01, dtype = np.float64)
  b = np.arange(-3, 3, 0.01, dtype = np.float64)

  xg, yg = np.meshgrid(a, b)

  z = (np.sin(xg) - yg)**2 + (xg - yg)**2

  xval = np.linspace(0, 5, 10)

  return x1, y1, x2, y2, xg, yg, z, xval  


def main():
  x = []	# координаты по х
  y = []	# координаты по у

  with open("coords.txt") as file:
  	for line in file:
  		xc, yc = line.split(" ")
  		x.append(dcm(xc))
  		y.append(dcm(yc))


  x1, y1, x2, y2, xg, yg, z, xval = test2()
  #x1, y1, x2, y2, xl, yl, xg, yg, z, xval = test4()
  #x1, y1, y2, y3, xg, yg, z, xval = test5()
  #x1, y1, x2, y2, xg, yg, z, xval = test6()

  fig, ax = plt.subplots()
  plt.grid()
  plt.gca().set_aspect('equal')
  ax.plot(x1, y1, label = 'F1', color = 'red', linewidth = 2)
  ax.plot(x2, y2, label = 'F2', color = 'yellow', linewidth = 2)
  #ax.plot(xl, yl, label = 'F3', color = 'blue', linewidth = 2)
  plt.legend(loc = 'upper left')

  residual = pl.contourf(xg, yg, z, levels = xval)
  plt.colorbar(residual, location = "left")

  plt.scatter(x, y, s = 40, color = "red")
  ax.plot(x, y)

  plt.show()

if __name__ == "__main__":
  main()