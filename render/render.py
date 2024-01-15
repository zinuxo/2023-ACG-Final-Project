import taichi as ti
import numpy as np
import time

ti.init(arch=ti.gpu)  # Try to run on GPU

res = (1080, 720)
window = ti.ui.Window("MPM 3D", res, vsync=True)

def render(x, colors):
    camera.track_user_inputs(window, movement_speed=0.03, hold_key=ti.ui.RMB)
    scene.set_camera(camera)

    scene.ambient_light((0, 0, 0))

    
    scene.particles(x, per_vertex_color = colors, radius=0.01)

    scene.point_light(pos=(0.5, 1.5, 0.5), color=(0.5, 0.5, 0.5))
    scene.point_light(pos=(0.5, 1.5, 1.5), color=(0.5, 0.5, 0.5))
    scene.point_light(pos=(0.5, -0.5, 0.5), color=(0.5, 0.5, 0.5))
    scene.point_light(pos=(0.5, -0.5, -0.5), color=(0.5, 0.5, 0.5))

    canvas.scene(scene)

canvas = window.get_canvas()
gui = window.get_gui()
scene = ti.ui.Scene()
camera = ti.ui.Camera()
# camera.position(-0.6,-1.4,-0.9)
# camera.position(0.7,1.4,2.5)
# camera.position(0.5,0.6,1.7)
# camera.position(0.5,0.5,-1)
# camera.position(0.5,-1,-1)
# camera.lookat(0.5, 0.0, 0.5)
# camera.position(0.5,1.8,-1.5)
# camera.position(0.9,1.3,1.5)
# camera.position(0.9,1,0.3)
# camera.position(0.5,2.5,0.5)
# camera.lookat(0.5, 0.2, 0.5)

# camera.position(1.5,1,-2)
# camera.position(3.5,1.8,1)
# camera.position(4.5,-2,1)
# camera.position(-1.5,2,1)
# camera.position(2.5,2,-2)
camera.position(0.5,4,1)
camera.lookat(1.5,1,1)
# camera.position(0.5,4,4)
camera.fov(45)
    
def main():

    # for frame_id in range(200,300):
    # for frame_id in range(450,500):
    # for frame_id in range(15,19):
    # for frame_id in range(70,90):
    # for frame_id in range(120,150):
    # for frame_id in range(100,160):
    # for frame_id in range(1450,1520):
    for frame_id in range(100):
        print("frame_id: ", frame_id)
        data = np.loadtxt(f'raw/{frame_id}.txt',dtype=np.float32,skiprows=1)
        L = len(data)
        x = ti.Vector.field(3, float, L)
        colors = ti.Vector.field(3, float, L)
        x.from_numpy(data[:, :3])
        colors.from_numpy(data[:, 3:])

        render(x, colors)
        window.show()
        # time.sleep(100)


if __name__ == '__main__':
    main()