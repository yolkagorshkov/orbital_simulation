import pygame
from math import sin, cos, sqrt, atan2, pi
from math import log as ln


class Vector:
    def __init__(self, a, b, flag='c'):
        if flag == 'c':
            self.x, self.y = a, b
        else:
            self.x, self.y = a * cos(b), a * sin(b)

    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y)

    def __iadd__(self, other):
        self.x += other.x
        self.y += other.y
        return self

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y)

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        return self

    def __mul__(self, other: float):
        return Vector(self.x * other, self.y * other)

    def __rmul__(self, other: float):
        return Vector(self.x * other, self.y * other)

    def __truediv__(self, other : float):
        if other == 0:
            return Vector(0, 0)
        return Vector(self.x / other, self.y / other)

    def __str__(self):
        return f'({self.x}; {self.y})'

    def __abs__(self) -> float:
        return sqrt(self.x ** 2 + self.y ** 2)

    def __neg__(self):
        return Vector(-self.x, -self.y)

    def __invert__(self): #нормированный вектор
        return self / abs(self)

    def angle(self) -> float:
        return atan2(self.x, self.y)

    def scalar_multiplication(self, other) -> float:
        return self.x * other.x + self.y * other.y

    def vector_multiplication(self, other) -> float:
        return self.x * other.y - self.y * other.x

    def angle_between_vectors(self, other):
        return atan2(self.vector_multiplication(other), self.scalar_multiplication(other))

    def rotate(self, theta):
        x_stroke = self.x * cos(theta) - self.y * sin(theta)
        y_stroke = self.x * cos(theta) - self.y * sin(theta)
        return Vector(x_stroke, y_stroke)

    @property
    def tup(self):
        return self.x, self.y


P0 = 100_000#Па
R = 8.31
M0 = .029
T0 = 273

C_d = .3
A = 60
m0 = 1160_000 #2487

first_stage_dry_mass = 96_000#74 + 22 на посадку
second_stage_full_mass = 347_000
second_stage_dry_mass = 25_000

first_stage_fuel = m0 - second_stage_full_mass - first_stage_dry_mass #2140

first_stage_engines_number = 12
second_stage_engines_number = 2 #РД-0120

first_stage_mode_one_consumption = 491.6 * first_stage_engines_number
first_stage_mode_two_consumption = 173.2 * first_stage_engines_number
second_stage_consumption = 420.35 * second_stage_engines_number

mu = 3.98e14
r0 = 6_378_000


def C_x(mach: float) -> float:
    if mach <= .65:
        return .2
    if .65 < mach <= 1.5:
        return 3 * mach / 8 - 17 / 400
    else:
        return 4 / 7 - 6 / 175 * mach


def mach_number(velocity: Vector, r: Vector) -> float:
    altitude = abs(r) - r0
    if altitude < 10_000:
        return abs(velocity) / (340 - altitude / 250)
    if 10_000 <= altitude <= 25_000:
        return abs(velocity) / 295
    else:
        return abs(velocity) / 300


def consumption(stage: int, mode : int, delta_t: float):
    if stage == 1:
        if mode == 1:
            return first_stage_mode_one_consumption * delta_t
        return first_stage_mode_two_consumption * delta_t
    return second_stage_consumption * delta_t


def p(r: Vector): #зависимость давления от высоты
    altitude = abs(r) - r0
    return P0 * 10 ** (-.06 * altitude / 1000)


def rho(r: Vector): #зависимость плотности от высоты
    return (p(r) * M0) / (R * T0)


def stage_thrust(r: Vector, stage: int, mode: int) -> float:
    if stage == 1:
        if mode == 1:
            return (2e6 + (6122 - p(r)) * pi * 1.2 ** 2) * first_stage_engines_number #степень расширения сопла 70
        return (7.845e5 + (519 - p(r)) * pi * 1.2 ** 2) * first_stage_engines_number #степень расширения сопла 170
    return (1.962e6 + (2961 - p(r)) * pi * 1.21 ** 2) * second_stage_engines_number #степень расширения сопла 86


def tsiolkovsky_delta_v(mass0, mass1, I):
    return I * ln(mass0 / mass1)


def gravity_acceleration(r: Vector) -> Vector:
    return -r * mu / abs(r) ** 3


def vertical_flight_thrust_acceleration(thrust, mass, r: Vector) -> Vector:
    return ~r * thrust / mass


def initial_pitch_over_maneuver_thrust_acceleration(thrust, mass, r: Vector) -> Vector:
    return (~r).rotate(pitch_angle) * thrust / mass


def gravity_turn_thrust_acceleration(thrust, mass, velocity: Vector) -> Vector:
    return ~velocity * thrust / mass


def drag_acceleration(velocity: Vector, r: Vector, mass) -> Vector:
    mach = mach_number(velocity, r)
    atmospheric_density = rho(r)
    return -(~velocity) * .5 * atmospheric_density * abs(velocity) ** 2 * A * C_x(mach) / mass


def thrust_total_acceleration(r: Vector, velocity: Vector, mass, thrust, delta_t, mode='stnd') -> Vector:
    global gravity_loss
    global drag_loss

    gravity_loss += gravity_acceleration(r) * delta_t * cos(velocity.angle_between_vectors(r))
    drag_loss += drag_acceleration(velocity, r, mass) * delta_t

    if mode == 'vert':
        return gravity_acceleration(r) + drag_acceleration(velocity, r, mass) + vertical_flight_thrust_acceleration(thrust, mass, r)
    if mode == 'piov':
        return gravity_acceleration(r) + drag_acceleration(velocity, r, mass) + initial_pitch_over_maneuver_thrust_acceleration(thrust, mass, r)
    else:
        return gravity_acceleration(r) + drag_acceleration(velocity, r, mass) + gravity_turn_thrust_acceleration(thrust, mass, velocity)


def no_thrust_total_acceleration(r: Vector, velocity: Vector, mass):
    return gravity_acceleration(r) + drag_acceleration(velocity, r, mass)


def thrust_making_step(r: Vector, velocity: Vector, step_length):
    global stage
    global first_stage_mode
    global mass
    global had_pitch_over
    global pitch
    if first_stage_mode == 1 and had_pitch_over:
        if stage_thrust(r, 1, 2) >= second_stage_trust_to_weight_ratio * abs(gravity_acceleration(r) * mass) / cos((velocity + starting_y_velocity).angle_between_vectors(r)):
            first_stage_mode = 2
            print('Переключение режима двигателей первой ступени!\n')
    if stage == 1:
        if mass <= m0 - first_stage_fuel:
            mass -= first_stage_dry_mass
            stage = 2
            print('Отделение первой ступени!\n')

    thrust = stage_thrust(r, stage, first_stage_mode)
    if not had_pitch_over:
        acceleration = thrust_total_acceleration(r, velocity, mass, thrust, step_length, 'vert')

    else:
        if pitch:
            acceleration = thrust_total_acceleration(r, velocity, mass, thrust, step_length, 'piov')
        else:
            acceleration = thrust_total_acceleration(r, velocity, mass, thrust, step_length, 'stnd')

    mass -= consumption(stage, first_stage_mode, step_length)

    delta_r = velocity * step_length + 0.5 * acceleration * step_length ** 2
    return (r + delta_r, velocity + acceleration * step_length)


def no_thrust_making_step(r: Vector, velocity: Vector, step_length):
    acceleration = no_thrust_total_acceleration(r, velocity, mass)
    delta_r = velocity * step_length + 0.5 * acceleration * step_length ** 2
    return (r + delta_r, velocity + acceleration * step_length)


latitude = 51.9 * pi / 180
starting_y_velocity = Vector(493 * cos(latitude), 0)

r = Vector(0, r0)
velocity = starting_y_velocity
stage = 1
first_stage_mode = 1
mass = m0
gravity_loss = Vector(0, 0)
drag_loss = Vector(0, 0)


sps = 100
printing_interval = 1
simulation_time = 660

pitch_time = 1.6
pitch_angle = -0.0344
pith_altitude = 160
second_stage_trust_to_weight_ratio = 1

had_pitch_over = False
pitch = False
pitch_begin = -1

coordinates = []

for i in range(simulation_time * sps):

    if mass <= 87_000 + second_stage_dry_mass:
        break

    if abs(velocity) * sin((-r).angle_between_vectors(velocity)) > sqrt(mu / abs(r)):
        break

    if first_stage_mode == 1:
        r, velocity = thrust_making_step(r, velocity - starting_y_velocity, 1 / sps)
        velocity += starting_y_velocity
    else:
        r, velocity = thrust_making_step(r, velocity, 1 / sps)

    if pitch and i >= pitch_begin + sps * pitch_time:
        pitch = False
    if not had_pitch_over and abs(r) - r0 >= pith_altitude:
        pitch_begin = i
        had_pitch_over = True
        pitch = True
    if i % (sps * printing_interval) == 0:
        if first_stage_mode == 1:
            coordinates.append((r, (0xff, 0, 0), abs(velocity)))
        elif stage == 1:
            coordinates.append((r, (0xff, 0xff, 0), abs(velocity)))
        else:
            coordinates.append((r, (0, 0xff, 0), abs(velocity)))

        '''print(f'Секунда: {i / sps}:\nv = {velocity}, alt = {(abs(r) - r0)}, mass = {mass}, normal ang = {r.angle() * 180 / pi}, velocity ang = {velocity.angle_between_vectors(r) * 180 / pi}')
        print(f'ang1 = {(r - Vector(0, r0)).angle_between_vectors(Vector(0, r0)) * 180 / pi}')'''
        print(f'orbital speed = {abs(velocity) * sin((-r).angle_between_vectors(velocity))}, delta alt speed = {abs(velocity) * cos((r).angle_between_vectors(velocity))}')
        '''if first_stage_mode == 1:
            print(f'fact velocity = {velocity - starting_y_velocity}, fact velocity angle = {(velocity - starting_y_velocity).angle_between_vectors(r) * 180 / pi}')
'''



print(gravity_loss, drag_loss)
print('\nНачало орбитального полета!\n')


orbital_simulation_time = 6000
orbital_sps = 100
orbital_printing_interval = 1

for i in range(orbital_simulation_time * orbital_sps):
    r, velocity = no_thrust_making_step(r, velocity, 1 / orbital_sps)

    if i % (orbital_sps * orbital_printing_interval) == 0:
        coordinates.append((r, (0xff, 0xff, 0xff), abs(velocity)))
        '''print(f'Секунда орбитального полета: {i / orbital_sps}:\nv = {velocity}, alt = {(abs(r) - r0)}, normal ang = {r.angle() * 180 / pi}')
        print(f'orbital speed = {abs(velocity) * sin((-r).angle_between_vectors(velocity))}, delta alt speed = {abs(velocity) * cos((r).angle_between_vectors(velocity))}')
'''

pygame.init()
screen = pygame.display.set_mode((1200, 1200))
pygame.display.set_caption("Orbital Flight Simulation")
clock = pygame.time.Clock()
running = True
font_arial = pygame.font.SysFont('arial', 36)

FPS = 60

scale = 6378
zero = Vector(0, 0)

motion = Vector(0, 0)

curr_t = 0
t_since_start = 0
time_warp = 60

while running:
    clock.tick(FPS)
    events = pygame.event.get()

    for event in events:
        if event.type == pygame.QUIT:
            running = False
        if event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 4:
                scale /= 1.016
            if event.button == 5:
                scale *= 1.016
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_d:
                motion.x = -r0 * 0.01 * (scale / 6378)
            if event.key == pygame.K_a:
                motion.x = r0 * 0.01 * (scale / 6378)
            if event.key == pygame.K_w:
                motion.y = -r0 * 0.01 * (scale / 6378)
            if event.key == pygame.K_s:
                motion.y = r0 * 0.01 * (scale / 6378)

        if event.type == pygame.KEYUP:
            if event.key == pygame.K_d or event.key == pygame.K_a:
                motion.x = 0
            if event.key == pygame.K_s or event.key == pygame.K_w:
                motion.y = 0

    zero += motion

    screen.fill((0, 0, 0))

    xc, yc = (zero * (1 / scale)).tup
    pygame.draw.circle(screen, (0x0, 0x0, 0x8a), (xc, 1200 - yc), r0 / scale)
    for i in range(min(len(coordinates) - 1, int(curr_t))):
        x1, y1 = ((coordinates[i][0] + zero) * (1 / scale)).tup
        x2, y2 = ((coordinates[i + 1][0] + zero) * (1 / scale)).tup
        pygame.draw.line(screen, coordinates[i][1], (x1, 1200 - y1), (x2, 1200 - y2), 5)

    curr_t += time_warp / FPS

    #pygame.draw.rect(screen, (0, 0, 0), ((0, 1160), (100, 1200)))
    try:
        alt = font_arial.render('altitude ' + str(round((abs(coordinates[int(curr_t)][0]) - 6378000)/1000, 1)), False, (0xff, 0xff, 0xff), None)
        vel = font_arial.render('velocity ' + str(round(abs(coordinates[int(curr_t)][2]), 1)), False, (0xff, 0xff, 0xff), None)
        screen.blit(alt, (0, 1160))
        screen.blit(vel, (0, 1120))
    except IndexError:
        pass
    pygame.display.flip()

pygame.quit()
