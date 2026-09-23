# import numpy as np
# import matplotlib.pyplot as plt

# # Parameters for planet and deformation
# radius = 1  # Radius of the planet
# deformation_factor = 0.2  # Factor by which the planet deforms due to tidal forces

# # Create a grid of angles (theta) to represent a circle (the planet)
# theta = np.linspace(0, 2 * np.pi, 100)
# x_circle = radius * np.cos(theta)
# y_circle = radius * np.sin(theta)

# # Apply tidal deformation: elongation along the x-axis and compression along y-axis
# x_deformed = (1 + deformation_factor) * x_circle  # Stretched along x-axis
# y_deformed = (1 - deformation_factor) * y_circle  # Compressed along y-axis

# # Plotting the original planet (circle) and the deformed planet
# fig, ax = plt.subplots(figsize=(6, 6))

# # Plot the original planet
# ax.plot(x_circle, y_circle, label="Original shape", linestyle='--', color='blue')

# # Plot the deformed planet due to tidal forces
# ax.plot(x_deformed, y_deformed, label="Deformed shape", color='red')

# # Add arrows to represent the deformation direction
# for angle in np.linspace(0, 2 * np.pi, 12):
#     dx = deformation_factor * np.cos(angle)
#     dy = -deformation_factor * np.sin(angle)
#     ax.arrow(np.cos(angle), np.sin(angle), dx, dy, head_width=0.05, color='green')

# # Add some labels and styling
# ax.set_aspect('equal')  # Ensure the aspect ratio is equal to maintain the circular shape
# ax.set_xlabel("X-axis (Tidal direction)")
# ax.set_ylabel("Y-axis (Perpendicular direction)")
# ax.set_title("Tidal Deformation of a Planet")
# ax.legend()

# plt.grid(True)
# plt.show()
####################################################################################
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Parameters for the planet and deformation
# radius = 1  # Radius of the planet
# deformation_factor_x = 0.3  # Stretch along x-axis (tidal axis)
# deformation_factor_y = 0.2  # Compression along y-axis (perpendicular to tidal axis)
# deformation_factor_z = 0.1  # Slight compression or stretch along the z-axis

# # Create spherical coordinates for the initial planet (theta: azimuthal, phi: polar)
# theta = np.linspace(0, 2 * np.pi, 100)
# phi = np.linspace(0, np.pi, 100)
# theta, phi = np.meshgrid(theta, phi)

# # Initial spherical coordinates (x, y, z) for the planet
# x_sphere = radius * np.sin(phi) * np.cos(theta)
# y_sphere = radius * np.sin(phi) * np.sin(theta)
# z_sphere = radius * np.cos(phi)

# # Apply tidal deformation
# x_deformed = (1 + deformation_factor_x) * x_sphere  # Stretch along x-axis
# y_deformed = (1 - deformation_factor_y) * y_sphere  # Compress along y-axis
# z_deformed = (1 - deformation_factor_z) * z_sphere  # Slight compression/stretch along z-axis

# # Plotting the original planet and the deformed planet
# fig = plt.figure(figsize=(10, 6))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the original planet (sphere)
# ax.plot_surface(x_sphere, y_sphere, z_sphere, color='lightblue', alpha=0.5, label='Original Planet')

# # Plot the deformed planet (ellipsoid due to tidal forces)
# ax.plot_surface(x_deformed, y_deformed, z_deformed, color='red', alpha=0.6)

# # Customizing the view
# ax.set_xlabel('X-axis (Tidal Axis)')
# ax.set_ylabel('Y-axis (Perpendicular Axis)')
# ax.set_zlabel('Z-axis (Vertical)')
# ax.set_title("3D Visualization of Tidal Deformation on a Planet")

# # Adjust the aspect ratio to make the plot look proportional
# ax.set_box_aspect([1, 1, 1])

# # Show plot
# plt.show()
##############################################################################
# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Parameters for the planet and deformation
# radius = 1  # Radius of the planet
# deformation_factor_y = 0.3  # Stretch along y-axis (tidal axis)
# deformation_factor_x = 0.2  # Compression along x-axis (perpendicular to tidal axis)
# deformation_factor_z = 0.1  # Slight compression or stretch along the z-axis

# # Create spherical coordinates for the initial planet (theta: azimuthal, phi: polar)
# theta = np.linspace(0, 2 * np.pi, 100)
# phi = np.linspace(0, np.pi, 100)
# theta, phi = np.meshgrid(theta, phi)

# # Initial spherical coordinates (x, y, z) for the planet
# x_sphere = radius * np.sin(phi) * np.cos(theta)
# y_sphere = radius * np.sin(phi) * np.sin(theta)
# z_sphere = radius * np.cos(phi)

# # Apply tidal deformation
# x_deformed = (1 - deformation_factor_x) * x_sphere  # Compress along x-axis
# y_deformed = (1 + deformation_factor_y) * y_sphere  # Stretch along y-axis
# z_deformed = (1 - deformation_factor_z) * z_sphere  # Slight compression/stretch along z-axis

# # Plotting the original planet and the deformed planet
# fig = plt.figure(figsize=(10, 6))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the original planet (sphere)
# ax.plot_surface(x_sphere, y_sphere, z_sphere, color='lightblue', alpha=0.5, label='Original Planet')

# # Plot the deformed planet (ellipsoid due to tidal forces along Y-axis)
# ax.plot_surface(x_deformed, y_deformed, z_deformed, color='red', alpha=0.6)

# # Adding centered x, y, z axes
# max_limit = radius * (1 + max(deformation_factor_x, deformation_factor_y, deformation_factor_z)) * 1.5
# ax.quiver(0, 0, 0, max_limit, 0, 0, color="black", arrow_length_ratio=0.05)  # X-axis
# ax.quiver(0, 0, 0, 0, max_limit, 0, color="black", arrow_length_ratio=0.05)  # Y-axis
# ax.quiver(0, 0, 0, 0, 0, max_limit, color="black", arrow_length_ratio=0.05)  # Z-axis

# # Customizing the view
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y-axis (Tidal Axis)')
# ax.set_zlabel('Z-axis')
# ax.set_title("3D Visualization of Tidal Deformation on a Planet (Deformation Along Y-axis)")

# # Adjust the aspect ratio to make the plot look proportional
# ax.set_box_aspect([1, 1, 1])

# # Set axis limits for better visualization
# ax.set_xlim([-max_limit, max_limit])
# ax.set_ylim([-max_limit, max_limit])
# ax.set_zlim([-max_limit, max_limit])

# # Show plot
# plt.show()
##############################################################################


# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Parameters for the planet and deformation
# radius = 2  # Radius of the planet
# deformation_factor_y = 0.3  # Stretch along y-axis (tidal axis)
# deformation_factor_x = 0.2  # Compression along x-axis (perpendicular to tidal axis)
# deformation_factor_z = 0.1  # Slight compression or stretch along the z-axis

# # Create spherical coordinates for the initial planet (theta: azimuthal, phi: polar)
# theta = np.linspace(0, 2 * np.pi, 100)
# phi = np.linspace(0, np.pi, 100)
# theta, phi = np.meshgrid(theta, phi)

# # Initial spherical coordinates (x, y, z) for the planet
# x_sphere = radius * np.sin(phi) * np.cos(theta)
# y_sphere = radius * np.sin(phi) * np.sin(theta)
# z_sphere = radius * np.cos(phi)

# # Apply tidal deformation
# x_deformed = (1 - deformation_factor_x) * x_sphere  # Compress along x-axis
# y_deformed = (1 + deformation_factor_y) * y_sphere  # Stretch along y-axis
# z_deformed = (1 - deformation_factor_z) * z_sphere  # Slight compression/stretch along z-axis

# # Plotting the original planet and the deformed planet
# fig = plt.figure(figsize=(25, 15))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the original planet (sphere) with reduced opacity
# ax.plot_surface(x_sphere, y_sphere, z_sphere, color='lightblue', alpha=0.4, label='Original Planet')

# # Plot the deformed planet (ellipsoid due to tidal forces along Y-axis) with reduced opacity
# ax.plot_surface(x_deformed, y_deformed, z_deformed, color='red', alpha=0.2)

# # Adding centered x, y, z axes
# max_limit = radius * (1 + max(deformation_factor_x, deformation_factor_y, deformation_factor_z)) * 1.5
# ax.quiver(0, 0, 0, max_limit, 0, 0, color="black", arrow_length_ratio=0.05, linewidth=2)  # X-axis
# ax.quiver(0, 0, 0, 0, max_limit, 0, color="black", arrow_length_ratio=0.05, linewidth=2)  # Y-axis
# ax.quiver(0, 0, 0, 0, 0, max_limit, color="black", arrow_length_ratio=0.05, linewidth=2)  # Z-axis

# # Plot deformation arrows with increased size and visibility
# n_arrows = 10  # Number of arrows to show deformation
# arrow_color = 'black'  # Color for the arrows
# arrow_size = 0.9  # Size of the arrows

# for i in range(n_arrows):
#     for j in range(n_arrows):
#         # Select points in spherical coordinates
#         x_start = x_sphere[i * 10, j * 10]
#         y_start = y_sphere[i * 10, j * 10]
#         z_start = z_sphere[i * 10, j * 10]
        
#         # Deformation vectors
#         dx = x_deformed[i * 10, j * 10] - x_start
#         dy = y_deformed[i * 10, j * 10] - y_start
#         dz = z_deformed[i * 10, j * 10] - z_start
        
#         # Plot arrows to show deformation direction with larger arrows and distinct color
#         ax.quiver(x_start, y_start, z_start, dx, dy, dz, color=arrow_color, 
#                   arrow_length_ratio=0.35, linewidth=2, length=arrow_size)

# # Customizing the view
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y-axis (Tidal Axis)')
# ax.set_zlabel('Z-axis')
# ax.set_title("3D Visualization of Tidal Deformation on a Planet (Deformation Along Y-axis)")

# # Adjust the aspect ratio to make the plot look proportional
# ax.set_box_aspect([1, 1, 1])

# # Set axis limits for better visualization
# ax.set_xlim([-max_limit, max_limit])
# ax.set_ylim([-max_limit, max_limit])
# ax.set_zlim([-max_limit, max_limit])

# # Show plot
# plt.show()


############################################################################## ROTATIONAL FLATT


# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Parameters for the planet and rotational flattening
# radius = 1  # Radius of the planet
# deformation_factor_z = 0.4  # Flattening along z-axis (polar flattening)
# deformation_factor_x = 0.1  # Slight expansion in the equatorial plane (x and y)

# # Create spherical coordinates for the initial planet (theta: azimuthal, phi: polar)
# theta = np.linspace(0, 2 * np.pi, 100)
# phi = np.linspace(0, np.pi, 100)
# theta, phi = np.meshgrid(theta, phi)

# # Initial spherical coordinates (x, y, z) for the planet
# x_sphere = radius * np.sin(phi) * np.cos(theta)
# y_sphere = radius * np.sin(phi) * np.sin(theta)
# z_sphere = radius * np.cos(phi)

# # Apply rotational flattening deformation
# x_deformed = (1 + deformation_factor_x) * x_sphere  # Expand along x-axis
# y_deformed = (1 + deformation_factor_x) * y_sphere  # Expand along y-axis
# z_deformed = (1 - deformation_factor_z) * z_sphere  # Flatten along z-axis

# # Plotting the original planet and the deformed planet
# fig = fig = plt.figure(figsize=(25, 15))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the original planet (sphere) with reduced opacity
# ax.plot_surface(x_sphere, y_sphere, z_sphere, color='lightblue', alpha=0.2, label='Original Planet')

# # Plot the deformed planet (ellipsoid due to rotational flattening) with reduced opacity
# ax.plot_surface(x_deformed, y_deformed, z_deformed, color='red', alpha=0.3)

# # Adding centered x, y, z axes
# max_limit = radius * (1 + max(deformation_factor_x, deformation_factor_z)) * 1.5
# ax.quiver(0, 0, 0, max_limit, 0, 0, color="black", arrow_length_ratio=0.05, linewidth=2)  # X-axis
# ax.quiver(0, 0, 0, 0, max_limit, 0, color="black", arrow_length_ratio=0.05, linewidth=2)  # Y-axis
# ax.quiver(0, 0, 0, 0, 0, max_limit, color="black", arrow_length_ratio=0.05, linewidth=2)  # Z-axis

# # Plot deformation arrows with increased size and visibility
# n_arrows = 10  # Number of arrows to show deformation
# arrow_color = 'black'  # Color for the arrows
# arrow_size = 0.8  # Size of the arrows

# for i in range(n_arrows):
#     for j in range(n_arrows):
#         # Select points in spherical coordinates
#         x_start = x_sphere[i * 10, j * 10]
#         y_start = y_sphere[i * 10, j * 10]
#         z_start = z_sphere[i * 10, j * 10]
        
#         # Deformation vectors
#         dx = x_deformed[i * 10, j * 10] - x_start
#         dy = y_deformed[i * 10, j * 10] - y_start
#         dz = z_deformed[i * 10, j * 10] - z_start
        
#         # Plot arrows to show deformation direction with larger arrows and distinct color
#         ax.quiver(x_start, y_start, z_start, dx, dy, dz, color=arrow_color, 
#                   arrow_length_ratio=0.35, linewidth=2, length=arrow_size)

# # Customizing the view
# ax.set_xlabel('X-axis (Equatorial Plane)')
# ax.set_ylabel('Y-axis (Equatorial Plane)')
# ax.set_zlabel('Z-axis (Polar Axis)')
# ax.set_title("3D Visualization of Rotational Flattening of a Planet")

# # Adjust the aspect ratio to make the plot look proportional
# ax.set_box_aspect([1, 1, 1])

# # Set axis limits for better visualization
# ax.set_xlim([-max_limit, max_limit])
# ax.set_ylim([-max_limit, max_limit])
# ax.set_zlim([-max_limit, max_limit])

# # Show plot
# plt.show()

##############################################################################



import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters for the planet and rotational flattening
radius = 1  # Radius of the planet
deformation_factor_z = 0.4  # Flattening along z-axis (polar flattening)
deformation_factor_x = 0.1  # Slight expansion in the equatorial plane (x and y)

# Create spherical coordinates for the initial planet (theta: azimuthal, phi: polar)
theta = np.linspace(0, 2 * np.pi, 100)
phi = np.linspace(0, np.pi, 100)
theta, phi = np.meshgrid(theta, phi)

# Initial spherical coordinates (x, y, z) for the planet
x_sphere = radius * np.sin(phi) * np.cos(theta)
y_sphere = radius * np.sin(phi) * np.sin(theta)
z_sphere = radius * np.cos(phi)

# Apply rotational flattening deformation
x_deformed = (1 + deformation_factor_x) * x_sphere  # Expand along x-axis
y_deformed = (1 + deformation_factor_x) * y_sphere  # Expand along y-axis
z_deformed = (1 - deformation_factor_z) * z_sphere  # Flatten along z-axis

# Plotting the original planet and the deformed planet
fig = fig = plt.figure(figsize=(25, 15))
ax = fig.add_subplot(111, projection='3d')

# Plot the original planet (sphere) with reduced opacity and latitude/longitude lines
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='lightblue', alpha=0.2, label='Original Planet')

# Add longitude and latitude lines to the original sphere
ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color='blue', alpha=0.3, linewidth=0.5)

# Plot the deformed planet (ellipsoid due to rotational flattening) with reduced opacity
ax.plot_surface(x_deformed, y_deformed, z_deformed, color='red', alpha=0.3)

# Add longitude and latitude lines to the deformed sphere
ax.plot_wireframe(x_deformed, y_deformed, z_deformed, color='red', alpha=0.5, linewidth=0.5)

# Adding centered x, y, z axes
max_limit = radius * (1 + max(deformation_factor_x, deformation_factor_z)) * 1.5
ax.quiver(0, 0, 0, max_limit, 0, 0, color="black", arrow_length_ratio=0.05, linewidth=2)  # X-axis
ax.quiver(0, 0, 0, 0, max_limit, 0, color="black", arrow_length_ratio=0.05, linewidth=2)  # Y-axis
ax.quiver(0, 0, 0, 0, 0, max_limit, color="black", arrow_length_ratio=0.05, linewidth=2)  # Z-axis

# Plot deformation arrows with increased size and visibility
n_arrows = 8  # Number of arrows to show deformation
arrow_color = 'black'  # Color for the arrows
arrow_size = 0.8  # Size of the arrows

for i in range(n_arrows):
    for j in range(n_arrows):
        # Select points in spherical coordinates
        x_start = x_sphere[i * 10, j * 10]
        y_start = y_sphere[i * 10, j * 10]
        z_start = z_sphere[i * 10, j * 10]
        
        # Deformation vectors
        dx = x_deformed[i * 10, j * 10] - x_start
        dy = y_deformed[i * 10, j * 10] - y_start
        dz = z_deformed[i * 10, j * 10] - z_start
        
        # Plot arrows to show deformation direction with larger arrows and distinct color
        ax.quiver(x_start, y_start, z_start, dx, dy, dz, color=arrow_color, 
                  arrow_length_ratio=0.35, linewidth=2, length=arrow_size)

# Customizing the view
ax.set_xlabel('X-axis (Equatorial Plane)')
ax.set_ylabel('Y-axis (Equatorial Plane)')
ax.set_zlabel('Z-axis (Polar Axis)')
ax.set_title("3D Visualization of Rotational Flattening with Latitude and Longitude Lines")

# Adjust the aspect ratio to make the plot look proportional
ax.set_box_aspect([1, 1, 1])

# Set axis limits for better visualization
ax.set_xlim([-max_limit, max_limit])
ax.set_ylim([-max_limit, max_limit])
ax.set_zlim([-max_limit, max_limit])

# Show plot
plt.show()

