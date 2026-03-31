# Unscented Kalman Filter (UKF) for Radar–LiDAR Sensor Fusion

This repository implements an **Unscented Kalman Filter (UKF)** in C++ for **nonlinear sensor fusion** using **RADAR** and **LiDAR** measurements.

The project is part of a **Sensor Fusion pipeline for autonomous driving**, where information from multiple sensors is combined to estimate the motion of surrounding vehicles more accurately than either sensor could achieve alone. Using a simulated highway environment, the system tracks multiple vehicles in real time and estimates their state under noisy measurements and nonlinear motion.

---

## 🎥 Demo – Highway Vehicle Tracking

![UKF Highway Tracking](ukf_highway_tracked.gif)

The demo shows a highway scenario with an ego vehicle and multiple surrounding traffic vehicles. The UKF fuses measurements from both sensors to track the position and motion of the observed vehicles over time.

### Visualization Details

- **Green vehicle**: ego vehicle
- **Blue vehicles**: surrounding traffic
- **Red markers**: LiDAR detections
- **Purple vectors**: RADAR measurements
- **Tracked state estimates**: UKF-based object state predictions
- **RMSE values**: estimation accuracy during runtime

---

## Project Overview

In autonomous driving, no single sensor is perfect:

- **LiDAR** provides accurate position measurements in Cartesian coordinates
- **RADAR** provides range, angle, and radial velocity, but in nonlinear polar coordinates

The challenge is to combine both measurement sources into a single reliable estimate of vehicle motion.

This project addresses that challenge by implementing an **Unscented Kalman Filter**, which is well suited for nonlinear systems because it propagates sigma points through the nonlinear process and measurement models instead of relying on linearization.

The implemented system:

- fuses **RADAR** and **LiDAR** data
- tracks **multiple vehicles**
- estimates the vehicle state in real time
- handles nonlinear motion and nonlinear measurements
- evaluates performance using **Root Mean Square Error (RMSE)**

---

## Sensor Fusion Context

This project is part of a broader **sensor fusion workflow** used in robotics and autonomous systems.

The UKF plays a critical role in that workflow by:

- combining multiple sensor streams into one consistent estimate
- reducing uncertainty from noisy measurements
- improving robustness when one sensor alone is insufficient
- enabling downstream tasks such as prediction, tracking, and decision-making

This makes the project directly relevant to:

- autonomous driving perception
- robotics tracking systems
- ADAS pipelines
- multi-sensor state estimation problems

---

## State Representation

The UKF estimates the following state variables for each tracked vehicle:

- **x-position**
- **y-position**
- **velocity**
- **yaw angle**
- **yaw rate**

This state representation allows the filter to model vehicle motion more realistically than a simple linear position-only tracker.

---

## Technical Approach

## 1. Initialization

The filter is initialized from the first incoming measurement. Depending on the sensor type:

- LiDAR measurements initialize position directly
- RADAR measurements are converted from polar to Cartesian coordinates before initialization

---

## 2. Prediction Step

At each time step, the UKF predicts the future vehicle state using a nonlinear motion model.

This includes:

- forming an **augmented state vector**
- generating **sigma points**
- propagating sigma points through the process model
- computing the predicted state mean
- computing the predicted covariance

This step allows the filter to estimate where the vehicle should be before receiving the next sensor measurement.

---

## 3. Update Step

When a new measurement arrives, the UKF corrects the predicted state.

### LiDAR Update
LiDAR measurements are linear in Cartesian space, so the update uses position information directly.

### RADAR Update
RADAR measurements are nonlinear and include:

- range
- bearing
- radial velocity

The UKF transforms predicted sigma points into RADAR measurement space and computes the correction without linearizing the model.

---

## 4. Multi-Vehicle Tracking

Each traffic vehicle in the highway simulation is assigned its own UKF instance.  
This allows the system to maintain separate state estimates for multiple moving objects simultaneously.

---

## Why UKF Instead of EKF?

The **Unscented Kalman Filter** is especially useful here because both the motion model and RADAR measurement model are nonlinear.

Compared with an Extended Kalman Filter (EKF), the UKF:

- avoids explicit Jacobian derivation
- handles stronger nonlinearities more gracefully
- often provides better estimation accuracy
- is more stable for nonlinear sensor fusion problems

---

## Results

The UKF successfully tracks surrounding vehicles in a noisy highway environment using fused RADAR and LiDAR measurements.

### Key observations

- position estimates remain stable over time
- vehicle motion is tracked smoothly through nonlinear maneuvers
- fusion improves tracking reliability compared with a single-sensor setup
- RMSE remains low enough to demonstrate effective estimation performance

This confirms that UKF-based sensor fusion is well suited for vehicle tracking in autonomous driving scenarios.

---

## What This Project Demonstrates

This project demonstrates practical understanding of:

- nonlinear state estimation
- sensor fusion
- RADAR measurement modeling
- LiDAR integration
- sigma point generation
- uncertainty propagation
- multi-object tracking
- performance evaluation using RMSE

---

## Tools and Environment

- **C++**
- **CMake**
- **Eigen**
- **PCL (Point Cloud Library)**
- **RADAR and LiDAR sensor simulation**

---

## Repository Structure

```text
unscented-kalman-filter-sensor-fusion/
│
├── src/
│   ├── ukf.cpp
│   ├── ukf.h
│   ├── tools.cpp
│   ├── tools.h
│   ├── main.cpp
│   ├── highway.h
│   └── measurement_package.h
│
├── ukf_highway_tracked.gif
├── CMakeLists.txt
├── README.md
└── LICENSE
