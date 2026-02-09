# Unscented Kalman Filter (UKF) – Radar & Lidar Sensor Fusion

This repository contains a **C++ implementation of an Unscented Kalman Filter (UKF)** for nonlinear **sensor fusion**, developed to track multiple vehicles in a simulated highway environment using noisy **RADAR** and **LIDAR** measurements.

The project demonstrates robust state estimation under nonlinear motion and measurement models, which are critical in **autonomous driving**, **robotics perception**, and **multi-target tracking systems**.

---

## 🎥 Demo – Highway Vehicle Tracking

![UKF Highway Tracking](ukf_highway_tracked.gif)


---

## 📌 Project Overview

- Full implementation of the **Unscented Kalman Filter**
- **Multi-vehicle tracking** in a three-lane highway scenario
- Sensor fusion using **RADAR** and **LIDAR**
- Estimated state vector:
  - Position (x, y)
  - Velocity
  - Yaw angle
  - Yaw rate
- Independent UKF instance for each tracked vehicle

---

## 🧠 Technical Approach

**Prediction**
- Augmented state and covariance generation
- Sigma point creation using the Unscented Transform
- Nonlinear motion model propagation

**Update**
- Separate update models for RADAR and LIDAR
- Kalman gain computation using cross-correlation matrices
- State and covariance correction

**Evaluation**
- Accuracy evaluated using **RMSE**
- Stable estimation during acceleration and lane-change maneuvers

---

## Build & Run


**Requirements**
- C++11 or later  
- CMake ≥ 3.5  
- Eigen  
- PCL (Point Cloud Library)

**Build**
```bash
mkdir build
cd build
cmake ..
make
./ukf_highway
---



## 🙏 Acknowledgement

This project was developed as part of the **Udacity Sensor Fusion Nanodegree**.  
The implementation, customization, and documentation reflect my own work and understanding of Unscented Kalman Filters and multi-sensor fusion.


