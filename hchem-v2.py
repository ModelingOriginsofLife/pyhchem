# Copyright (C) 2014 nineties
# $Id: sim.py 2014-07-23 21:39:31 nineties $

#= A clone of Tim Hutton's artificial chemistry simulator. =

import numpy as np
import numpy.linalg as la
import pygame, inputbox, pickle
import re, sys, time, math

class Rule:
    def __init__(self, filename):
        self.cnt = 0
        self.num = None
        self.fill = []
        self.types = []
        self.map = {}
        self.color_count = 0
        self.colors = []
        self.colormap = {}
        self.wildcards = []
        self.name = []
        self.ruleb = {} # Rules for bounded pair
        self.ruleu = {} # Rules for unbounded pair
        self.rule_texts = []
        self.colortable = [
            (0, 0, 0), (255, 0, 0), (0, 255, 0), (0, 0, 255),
            (255, 255, 0), (255, 0, 255), (0, 255, 255),
            (255, 127, 0), (255, 0, 127), (127, 255, 0), (127, 0, 255),
            (0, 255, 127), (0, 127, 255)
        ]
        self.parse(filename)

    def gen_color(self, a):
        self.color_count += 1
        return self.colortable[(self.color_count - 1) % len(self.colortable)]

    def to_index(self, t, n):
        term = t + n
        if not term in self.map:
            self.map[term] = self.cnt
            self.colors.append(self.colormap[t])
            self.name.append(term)
            self.cnt += 1
        return self.map[term]

    def get_index(self, term):
        return self.map[term]

    def get_name(self, idx):
        return self.name[idx]

    def parse_expr(self, str):
        if "-" in str:
            bnd = True
            M0, M1 = str.split("-")
        else:
            bnd = False
            M0, M1 = str.split()
        p = re.search(r'([a-zA-Z]+)(\d+)', M0)
        q = re.search(r'([a-zA-Z]+)(\d+)', M1)
        return (p.group(1), p.group(2), q.group(1), q.group(2), bnd)

    def add_rule(self, L0, l0, L1, l1, lbnd, R0, r0, R1, r1, rbnd, p):
        LL0 = self.to_index(L0, l0)
        LL1 = self.to_index(L1, l1)
        RR0 = self.to_index(R0, r0)
        RR1 = self.to_index(R1, r1)
        if lbnd:
            if not (LL0, LL1) in self.ruleb:
                self.ruleb[(LL0, LL1)] = []
                if LL0 != LL1:
                    self.ruleb[(LL1, LL0)] = []
            self.ruleb[(LL0, LL1)].append((RR0, RR1, rbnd, p))
            if LL0 != LL1:
                self.ruleb[(LL1, LL0)].append((RR1, RR0, rbnd, p))
        else:
            if not (LL0, LL1) in self.ruleu:
                self.ruleu[(LL0, LL1)] = []
                if LL0 != LL1:
                    self.ruleu[(LL1, LL0)] = []
            self.ruleu[(LL0, LL1)].append((RR0, RR1, rbnd, p))
            if LL0 != LL1:
                self.ruleu[(LL1, LL0)].append((RR1, RR0, rbnd, p))

    def parse_rule(self, line):
        try:
            lhs, rhs = line.split("->")
        except:
            return
        self.rule_texts.append(line)
        prob = 1.0
        if ":" in rhs:
            rhs, p = rhs.split(":")
            prob = eval(p.strip())
        L0, l0, L1, l1, lbnd = self.parse_expr(lhs.strip())
        R0, r0, R1, r1, rbnd = self.parse_expr(rhs.strip())
        if L0 in self.wildcards and L1 in self.wildcards:
            if L0 == L1:
                for t in self.types:
                    _L0 = t
                    _L1 = t
                    if R0 == L0: _R0 = t
                    else:        _R0 = R0
                    if R1 == L0: _R1 = t
                    else:        _R1 = R1
                    self.add_rule(_L0, l0, _L1, l1, lbnd, _R0, r0, _R1, r1, rbnd, prob)
            else:
                for t0 in self.types:
                    for t1 in self.types:
                        if l0 == l1 and t0 > t1: continue
                        _L0 = t0
                        _L1 = t1
                        if R0 == L0:   _R0 = t0
                        elif R0 == L1: _R0 = t1
                        else:          _R0 = R0
                        if R1 == L0:   _R1 = t0
                        elif R1 == L1: _R1 = t1
                        else:          _R1 = R1
                        self.add_rule(_L0, l0, _L1, l1, lbnd, _R0, r0, _R1, r1, rbnd, prob)
        elif L0 in self.wildcards or L1 in self.wildcards:
            if L0 in self.wildcards:
                for t in self.types:
                    _L0 = t
                    if R0 == L0: _R0 = t
                    else:        _R0 = R0
                    if R1 == L0: _R1 = t
                    else:        _R1 = R1
                    self.add_rule(_L0, l0, L1, l1, lbnd, _R0, r0, _R1, r1, rbnd, prob)
            else:
                for t in self.types:
                    _L1 = t
                    if R0 == L1: _R0 = t
                    else:        _R0 = R0
                    if R1 == L1: _R1 = t
                    else:        _R1 = R1
                    self.add_rule(L0, l0, _L1, l1, lbnd, _R0, r0, _R1, r1, rbnd, prob)
        else:
            self.add_rule(L0, l0, L1, l1, lbnd, R0, r0, R1, r1, rbnd, prob)

    def add_type(self, t):
        self.types.append(t)
        self.colormap[t] = self.gen_color(t)

    def setup_types(self, str):
        lhs, rhs = str.split(":")
        for t in rhs.split(","):
            self.add_type(t.strip())

    def setup_wildcards(self, str):
        lhs, rhs = str.split(":")
        for w in rhs.split(","):
            self.wildcards.append(w.strip())

    def setup_fill(self, str):
        lhs, rhs = str.split(":")
        for decl in rhs.split(","):
            t, p = decl.strip().split(" ")
            self.fill.append((t, eval(p)))

    def parse(self, filename):
        f = open(filename, "r")
        while True:
            line = f.readline()
            if not line: break
            if line[0] == '#': continue
            line = line.strip()
            if line.find("type") == 0:
                self.setup_types(line)
            elif line.find("wildcard") == 0:
                self.setup_wildcards(line)
            elif line.find("number of particles") == 0:
                self.num = int(line.split(":")[1].strip())
            elif line.find("fill") == 0:
                self.setup_fill(line)
            elif "->" in line:
                self.parse_rule(line)
        f.close()

    def check(self, L0, L1, bound):
        if bound:
            if not (L0, L1) in self.ruleb:
                return None
            return self.ruleb[(L0, L1)]
        else:
            if not (L0, L1) in self.ruleu:
                return None
            return self.ruleu[(L0, L1)]

class HChem:
    # n   : number of particles
    # r   : radious of particles
    # v0  : initial velocity of particles
    # dt  : duration of one frame
    # k   : strength of bounds
    # w,h : width and height of the universe
    # seed: random seed
    def __init__(self, filename, n = 1000, r = 10, v0 = None, dt = 0.1,
                 width = 1200, height = 700, bucket_size = None, seed=None):
        self.rule = Rule(filename)
        if seed: np.random.seed(seed)
        if v0 == None: v0 = r
        if bucket_size == None: bucket_size = r

        self.n = n
        if self.rule.num: self.n = self.rule.num
        self.r = r
        self.dt = dt
        self.w = width
        self.h = height

        # Initialize positions of particles
        self.pos = np.zeros((self.n, 2))
        self.pos[:,0] = np.random.uniform(r, width-r, self.n)
        self.pos[:,1] = np.random.uniform(r, height-r, self.n)

        # Initialize velocities of particles
        direction = np.random.uniform(0, 2*np.pi, self.n)
        self.vel = np.zeros((self.n, 2))
        self.vel[:,0] = v0*np.cos(direction)
        self.vel[:,1] = v0*np.sin(direction)

        # Initialize types
        self.types = np.zeros(self.n, dtype=int)
        for k in xrange(self.n):
            p = np.random.uniform(0, 1)
            q = 0
            for (t, r) in self.rule.fill:
                q += r
                if p < q:
                    self.types[k] = self.rule.get_index(t)
                    break
        
        for k in xrange(self.n):
            if self.pos[k, 0]  > self.w/2:
                self.types[k] = self.rule.get_index("a1")
        for k in xrange(35):
            self.pos[k, 0] = self.w/2
            self.pos[k, 1] = k*2*self.r + self.r
            self.types[k] = self.rule.get_index("m0")

        self.fixed = np.zeros(self.n, dtype=int)

        # self.bounds[i] == list of indexes of particles which is bound to i.
        self.bounds = np.zeros(self.n, dtype=object)
        for i in xrange(self.n): self.bounds[i] = []

        # Initialize buckets for compute_collision detection
        self.bucket_size = bucket_size
        self.nbx = int(math.ceil(float(width)/bucket_size))
        self.nby = int(math.ceil(float(height)/bucket_size))
        self.buckets = np.zeros((self.nbx, self.nby), dtype=object)

    def load_rule(self, fname):
        self.rule = Rule(fname)

    def bucket_index(self, x):
        return (min(max(int(x[0]/self.bucket_size), 0), self.nbx-1),
                min(max(int(x[1]/self.bucket_size), 0), self.nby-1))

    def init_bucket(self):
        for i in xrange(self.nbx):
            for j in xrange(self.nby):
                self.buckets[i,j] = []
        for k in xrange(self.n):
            i, j = self.bucket_index(self.pos[k, :])
            self.buckets[i, j].append(k)

    def add_impulse_from_walls(self):
        r = self.r
        for k in xrange(self.n):
            x = self.pos[k, 0]
            y = self.pos[k, 1]
            vx = self.vel[k, 0]
            vy = self.vel[k, 1]
            if (x < r and vx < 0) or (x > self.w-r and vx > 0):
                self.vel[k, 0] += -2*self.vel[k, 0]
                self.vel[k, 1] += 0
            if (y < r and vy < 0) or (y > self.h-r and vy > 0):
                self.vel[k, 0] += 0
                self.vel[k, 1] += -2*self.vel[k, 1]

    def update_state_of_particle_pair(self, k, l):
        if l in self.bounds[k]:
            # bound pair
            rules = self.rule.check(self.types[k], self.types[l], True)
            if rules:
                for r in rules:
                    p = r[3]
                    if np.random.uniform(0, 1) < p:
                        self.types[k] = r[0]
                        self.types[l] = r[1]
                        if not r[2]:
                            self.bounds[k].remove(l)
                            self.bounds[l].remove(k)
                            return False
                        return True
            return True
        else:
            # unbound pair
            rules = self.rule.check(self.types[k], self.types[l], False)
            if rules:
                for r in rules:
                    p = r[3]
                    if np.random.uniform(0, 1) < p:
                        self.types[k] = r[0]
                        self.types[l] = r[1]
                        if r[2]:
                            self.bounds[k].append(l)
                            self.bounds[l].append(k)
                            return True
                        return False
            return False

    def add_impulse_between_unbound_pair(self, k, l):
        rx = self.pos[k,:] - self.pos[l,:]
        d2 = np.sum(rx*rx)
        if d2 > 4*self.r*self.r: return
        rv = self.vel[k,:] - self.vel[l,:]
        c  = rx.dot(rv)
        if c >= 0: return
        if self.update_state_of_particle_pair(k, l):
            return
        d = math.sqrt(d2)
        fdt = -c*rx/d2
        self.vel[k,:] += fdt
        self.vel[l,:] -= fdt

    def add_impulse_between_bound_pair(self, k, l, rx, rv, d2):
        d = math.sqrt(d2)
        n = rx/d
        c = rx.dot(rv)
        #ldt = -(2*c + 2*(d2-4*self.r*self.r))/(8*d2)
        #self.vel[k,:] += 2*ldt*rx
        #self.vel[l,:] -= 2*ldt*rx
        if (d < 2*self.r and c < 0) or (d > 3*self.r and c > 0):
            ldt = -n.dot(rv)
            self.vel[k,:] += ldt*n
            self.vel[l,:] -= ldt*n

    def add_impulse_between_particles_sub(self, k, i, j):
        if i < 0 or j < 0 or i >= self.nbx or j >= self.nby: return
        for l in self.buckets[i, j]:
            if k >= l: continue
            self.add_impulse_between_unbound_pair(k, l)

    def add_impulse_between_particles(self):
        r = self.r

        # add impulse between unbound pairs
        for k in xrange(self.n):
            i,j = self.bucket_index(self.pos[k, :])
            self.add_impulse_between_particles_sub(k, i-1, j)
            self.add_impulse_between_particles_sub(k, i-1, j-1)
            self.add_impulse_between_particles_sub(k, i-1 ,j+1)
            self.add_impulse_between_particles_sub(k, i, j+1)
            self.add_impulse_between_particles_sub(k, i+1, j-1)
            self.add_impulse_between_particles_sub(k, i+1, j)
            self.add_impulse_between_particles_sub(k, i+1, j+1)

    def add_impulse_between_bound_particles(self):
        # add impulse between bound pairs
        for k in xrange(self.n):
            for l in self.bounds[k]:
                if k >= l: continue
                rx = self.pos[k,:] - self.pos[l,:]
                rv = self.vel[k,:] - self.vel[l,:]
                d2 = np.sum(rx*rx)
                self.add_impulse_between_bound_pair(k, l, rx, rv, d2)

    def add_impulse_from_bond(self, k, l, p, d, e):
        rx = self.pos[p,:] - self.pos[l,:]
        c  = rx.dot(e)
        if c < 0 or c > d: return
        h = rx - c*e
        if np.sum(h*h) > 1.1*self.r*self.r:
            return
        if self.vel[p,:].dot(h) > 0:
            return
        v  = self.vel[p, :]
        vh = v - v.dot(e)*e
        r  = c/d
        self.vel[p, :] += -2*vh
        self.vel[k, :] += r*2*vh
        self.vel[l, :] += (1-r)*2*vh

    def add_impulse_from_bonds_sub(self, k, l, i, j, d, e):
        for p in self.buckets[i, j]:
            if k == p or l == p: continue
            self.add_impulse_from_bond(k, l, p, d, e)

    def add_impulse_from_bonds(self):
        for k in xrange(self.n):
            for l in self.bounds[k]:
                if k >= l: continue
                rx = self.pos[k, :] - self.pos[l, :]
                d  = math.sqrt(np.sum(rx*rx))
                e  = rx/d
                i1, j1 = self.bucket_index(self.pos[k, :])
                i2, j2 = self.bucket_index(self.pos[l, :])
                if i1 > i2: i1, i2 = i2, i1
                if j1 > j2: j1, j2 = j2, j1
                i1 = max(i1-1, 0)
                i2 = min(i2+1, self.nbx-1)
                j1 = max(j1-1, 0)
                j2 = min(j2+1, self.nby-1)
                for i in xrange(i1, i2+1):
                    for j in xrange(j1, j2+1):
                        self.add_impulse_from_bonds_sub(k, l, i, j, d, e)


    def update(self):
        self.init_bucket()
        self.add_impulse_from_bonds()
        self.add_impulse_from_walls()
        self.add_impulse_between_bound_particles()
        self.add_impulse_between_particles()

        for k in xrange(self.n):
            if self.fixed[k] == 1:
                self.vel[k,:] = 0
        self.pos += self.vel*self.dt
        
    def total_energy(self):
        return np.sum(self.vel*self.vel)

    def save(self, fname):
        with open(fname, "w") as f:
            pickle.dump(self, f)

    @classmethod
    def load(self, fname):
        try:
            with open(fname, "r") as f:
                return pickle.load(f)
        except:
            pass

class HChemViewer:
    RED   = (255, 0, 0)
    BLUE  = (0, 0, 255)
    WHITE = (255, 255, 255)
    BLACK = (0, 0, 0)
    INFO = [
    "(P) play/pause, (F) stepwise, (Q) quit, (S) save snapshot, (L) load snapshot, (R) load rule, (T) show/hide particle types, (up) speed up, (down) slow down",
    "(left drag) move particle, (right drag) bind particles, (double click) change particle type"
    ]

    def __init__(self, sim, w = None, h = None):
        if w == None: w = sim.w
        if h == None: h = sim.h

        self.sim = sim
        pygame.init()
        self.screen = pygame.display.set_mode((w, h)
            #, pygame.DOUBLEBUF | pygame.FULLSCREEN | pygame.HWSURFACE
            )
        pygame.display.set_caption("Tim Hutton's Artifical Chemistry")
        self.fontsize = 18
        self.font = pygame.font.SysFont(None, self.fontsize)
        info_texts = self.INFO + sim.rule.rule_texts
        self.info = map(lambda text: self.font.render(text, False, self.BLUE),
                info_texts)

        self.speed = 10

        # For events
        self.play          = False
        self.stepwise      = False
        self.dragged       = None
        self.move          = False
        self.bound         = False
        self.display_types = False
        self.prev_lclick   = time.time()

    def get_clicked(self):
        for k in xrange(self.sim.n):
            d2 = np.sum((self.sim.pos[k,:] - pygame.mouse.get_pos())**2)
            if d2 < self.sim.r**2:
                return k
                break
        return None

    def check_event(self):
        for event in pygame.event.get():
            if not self.dragged and event.type == pygame.KEYDOWN:
                key = pygame.key.get_pressed()
                if key[pygame.K_q]:
                    sys.exit()
                if key[pygame.K_p]:
                    self.play = not self.play
                if key[pygame.K_s]:
                    fname = inputbox.ask(self.screen, "filename")
                    if fname: self.sim.save(fname)
                if key[pygame.K_l]:
                    fname = inputbox.ask(self.screen, "filename")
                    if fname: self.sim = HChem.load(fname)
                if key[pygame.K_t]:
                    self.display_types = not self.display_types
                if key[pygame.K_f]:
                    self.stepwise = True
                    self.play = True
                if key[pygame.K_r]:
                    fname = inputbox.ask(self.screen, "filename")
                    self.sim.load_rule(fname)
                    info_texts = self.INFO + sim.rule.rule_texts
                    self.info = map(lambda text: self.font.render(text, False,
                        self.BLUE), info_texts)
                    self.play = False
            elif not self.dragged and event.type == pygame.MOUSEBUTTONDOWN:
                self.play = False
                clicked = self.get_clicked()
                l,m,r = pygame.mouse.get_pressed()

                # Detect double click
                t = time.time()
                double_click = False
                if t - self.prev_lclick < 1.0/3:
                    double_click = True
                self.prev_lclick = t

                if l and double_click:
                    t = inputbox.ask(self.screen, "type")
                    try:
                        self.sim.types[clicked] = self.sim.rule.get_index(t)
                    except:
                        pass
                elif r and double_click:
                    self.sim.fixed[clicked] = not self.sim.fixed[clicked]
                elif clicked:
                    self.dragged = clicked
                    if l: self.move = True
                    elif r: self.bound = True
            elif self.dragged and event.type == pygame.MOUSEMOTION:
                l,m,r = pygame.mouse.get_pressed()
                if l:
                    self.sim.pos[self.dragged,:] = pygame.mouse.get_pos()
            elif self.dragged and event.type == pygame.MOUSEBUTTONUP:
                if self.bound:
                     clicked = self.get_clicked()
                     if clicked and self.dragged != clicked and\
                             not (clicked in self.sim.bounds[self.dragged]):
                         self.sim.bounds[self.dragged].append(clicked)
                         self.sim.bounds[clicked].append(self.dragged)
                self.move    = False
                self.bound   = False
                self.dragged = None
            elif event.type == pygame.QUIT:
                sys.exit()

    def loop(self):
        iteration = 0
        screen = self.screen
        while True:
            sim = self.sim
            n   = sim.n
            r   = sim.r
            if self.play:
                iteration += 1
                sim.update()

            if self.stepwise:
                self.play = False
                self.stepwise = False

            pos = sim.pos

            screen.fill(self.WHITE)
            # Draw particles
            for k in xrange(n):
                pygame.draw.circle(screen, sim.rule.colors[sim.types[k]],
                        (int(pos[k,0]), int(pos[k,1])), r, 1)

            if self.display_types:
                for k in xrange(sim.n):
                    t = sim.rule.get_name(sim.types[k])
                    text = self.font.render(t, False, self.BLACK)
                    rect = text.get_rect()
                    rect.centerx = pos[k,0]
                    rect.centery = pos[k,1]
                    self.screen.blit(text, rect)

            # Draw bounds
            for k in xrange(n):
                for l in sim.bounds[k]:
                    if k >= l: continue
                    pygame.draw.line(screen, self.BLACK, pos[k,:], pos[l,:])

            # Other info
            if self.bound:
                pygame.draw.line(screen, self.BLACK,
                    pos[self.dragged,:], pygame.mouse.get_pos())

            y = 10
            for i in self.info:
                self.screen.blit(i, (10, y))
                y += self.fontsize
            text = self.font.render(
                    "time = " + str(iteration*sim.dt),
                    False, self.BLUE)
            self.screen.blit(text, (10, y))
            energy = sim.total_energy()
            text = self.font.render(
                    "energy = " + str(energy),
                    False, self.BLUE)
            self.screen.blit(text, (10, y + self.fontsize))
            self.check_event()
            pygame.display.update()

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        sim = HChem(sys.argv[1])
    else:
        sim = HChem("membrane.txt")
    HChemViewer(sim).loop()
