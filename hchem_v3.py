# -*- coding: utf-8 -*-
# Copyright (C) 2014 nineties
# $Id: pyhchem_v3.py 2014-08-05 15:19:41 nineties $

#= A clone of Tim Hutton's artificial chemistry simulator. =

import numpy as np
import Tkinter as tk
import tkFileDialog
import pygame
import math, re, time, sys, os

class HChemSimulator:
    def __init__(self, n, types, init, rules, wildcards = [],
            width = 1300, height = 650, radious = 10, v0 = 10, dt = 0.25):
        self.n         = n
        self.types     = types
        self.wildcards = wildcards
        self.init      = init
        self.rules     = rules
        self.ruleu     = {} # rules for unbound pair
        self.ruleb     = {} # rules for bound pair
        self.w         = width
        self.h         = height
        self.r         = radious
        self.l         = 3*self.r
        self.v0        = v0
        self.dt        = dt
        self.setup_particles()
        # Warnup
        for k in range(5):
            self.update()
        self.setup_rules()
        print self.ruleu
        print self.ruleb

    # Convert "a0" to "a" and 0
    def parse_type(self, t):
        p = re.search(r'([a-zA-Z]+)(\d+)', t)
        return p.group(1), int(p.group(2))

    def setup_particles(self):
        # Initialize time
        self.t = 0

        # Initialize Mass
        self.mass = np.ones(self.n)

        # Initialize Positions
        self.pos = np.zeros((self.n, 2))
        self.pos[:,0] = np.random.uniform(self.r, self.w-self.r, self.n)
        self.pos[:,1] = np.random.uniform(self.r, self.h-self.r, self.n)

        # Initialize Velocities
        self.vel = np.zeros((self.n, 2))
        direction = np.random.uniform(0, 2*np.pi, self.n)
        self.vel[:,0] = self.v0 * np.cos(direction)
        self.vel[:,1] = self.v0 * np.sin(direction)

        # Initialize Bonds
        self.bonds = np.zeros(self.n, dtype=object)
        for k in xrange(self.n): self.bonds[k] = []

        # Initialize types and states
        self.typ   = np.zeros(self.n, dtype = object)
        self.state = np.zeros(self.n, dtype = int)

        for k in xrange(self.n):
            p = np.random.uniform(0, 1)
            q = 0
            for (t, r) in self.init:
                q += r
                if p <= q:
                    t, s = self.parse_type(t)
                    self.typ[k] = t
                    self.state[k] = s
                    break

        # Array of visited flags for DFS
        self.visited = np.zeros(self.n, dtype = int)

        # Initialize buckets for compute_collision detection
        self.bucket_size = 2*self.r
        self.nbx = int(math.ceil(float(self.w)/self.bucket_size))
        self.nby = int(math.ceil(float(self.h)/self.bucket_size))
        self.buckets = np.zeros((self.nbx, self.nby), dtype=object)

    def setup_rules(self):
        for rule in self.rules:
            if type(rule) != tuple:
                rule = (rule, 1.0)
            self.add_new_rule(rule[0], rule[1])

    def add_rule(self,lt0,ls0,lt1,ls1,lbnd,rt0,rs0,rt1,rs1,rbnd,p):
        if lbnd:
            if not (lt0,ls0,lt1,ls1) in self.ruleb:
                self.ruleb[(lt0,ls0,lt1,ls1)] = []
                if lt0 != lt1 and ls0 != ls1:
                    self.ruleb[(lt1,ls1,lt0,ls0)] = []
            self.ruleb[(lt0,ls0,lt1,ls1)].append((rt0,rs0,rt1,rs1,rbnd,p))
            if lt0 != lt1 and ls0 != ls1:
                self.ruleb[(lt1,ls1,lt0,ls0)].append((rt0,rs0,rt1,rs1,rbnd,p))
        else:
            if not (lt0,ls0,lt1,ls1) in self.ruleu:
                self.ruleu[(lt0,ls0,lt1,ls1)] = []
                if lt0 != lt1 and ls0 != ls1:
                    self.ruleu[(lt1,ls1,lt0,ls0)] = []
            self.ruleu[(lt0,ls0,lt1,ls1)].append((rt0,rs0,rt1,rs1,rbnd,p))
            if lt0 != lt1 and ls0 != ls1:
                self.ruleu[(lt1,ls1,lt0,ls0)].append((rt1,rs1,rt0,rs0,rbnd,p))

    def add_new_rule(self, rule, prob):
        m = re.search(r'([a-zA-Z]+)(\d+)(-|\s+)([a-zA-Z]+)(\d+)\s*->\s*([a-zA-Z]+)(\d+)(-|\s+)([a-zA-Z]+)(\d+)', rule)
        lt0  = m.group(1)
        ls0  = int(m.group(2))
        lbnd = m.group(3) == '-'
        lt1  = m.group(4)
        ls1  = int(m.group(5))
        rt0  = m.group(6)
        rs0  = int(m.group(7))
        rbnd = m.group(8) == '-'
        rt1  = m.group(9)
        rs1  = int(m.group(10))
        self.add_rule(lt0, ls0, lt1, ls1, lbnd, rt0, rs0, rt1, rs1, rbnd, prob)

    def check(self, k, l, bound):
        lt0 = self.typ[k]
        ls0 = self.state[k]
        lt1 = self.typ[l]
        ls1 = self.state[l]
        if bound:
            if not (lt0, ls0, lt1, ls1) in self.ruleb:
                return None
            return self.ruleb[(lt0, ls0, lt1, ls1)]
        else:
            if not (lt0, ls0, lt1, ls1) in self.ruleu:
                return None
            return self.ruleu[(lt0, ls0, lt1, ls1)]

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
        if l in self.bonds[k]:
            # bound pair
            rules = self.check(k, l, True)
            if rules:
                for r in rules:
                    p = r[5]
                    if np.random.uniform(0, 1) < p:
                        self.typ[k] = r[0]
                        self.state[k] = r[1]
                        self.typ[l] = r[2]
                        self.state[l] = r[3]
                        if not r[4]:
                            self.bonds[k].remove(l)
                            self.bonds[l].remove(k)
                            return False
                        return True
            return True
        else:
            # unbound pair
            rules = self.check(k, l, False)
            if rules:
                for r in rules:
                    p = r[5]
                    if np.random.uniform(0, 1) < p:
                        self.typ[k] = r[0]
                        self.state[k] = r[1]
                        self.typ[l] = r[2]
                        self.state[l] = r[3]
                        if r[4]:
                            self.bonds[k].append(l)
                            self.bonds[l].append(k)
                            return True
                        return False
            return False

    def add_impulse_between_unbound_pair(self, k, l, rx, rv, d2):
        if self.update_state_of_particle_pair(k, l):
            return
        d = math.sqrt(d2)
        n = rx/d
        ldt = -n.dot(rv)
        self.vel[k,:] += ldt*n
        self.vel[l,:] -= ldt*n

    def add_impulse_between_bound_pair(self, k, l, rx, rv, d2):
        d = math.sqrt(d2)
        n = rx/d
        c = rx.dot(rv)
        if (d < 2*self.r and c < 0) or (d > self.l and c > 0):
            ldt = -n.dot(rv)
            self.vel[k,:] += ldt*n
            self.vel[l,:] -= ldt*n

    def add_impulse_between_particles_sub(self, k, i, j):
        if i < 0 or j < 0 or i >= self.nbx or j >= self.nby: return
        for l in self.buckets[i, j]:
            if k >= l: continue
            rx = self.pos[k,:] - self.pos[l,:]
            rv = self.vel[k,:] - self.vel[l,:]
            if rx.dot(rv) >= 0: continue
            d2 = np.sum(rx*rx)
            if d2 > 4*self.r*self.r: continue
            self.add_impulse_between_unbound_pair(k, l, rx, rv, d2)

    def add_impulse_between_particles(self):
        r = self.r

        # add impulse between unbound pairs
        for k in xrange(self.n):
            i,j = self.bucket_index(self.pos[k, :])
            self.add_impulse_between_particles_sub(k, i-1, j)
            self.add_impulse_between_particles_sub(k, i-1, j-1)
            self.add_impulse_between_particles_sub(k, i-1 ,j+1)
            self.add_impulse_between_particles_sub(k, i, j-1)
            self.add_impulse_between_particles_sub(k, i, j)
            self.add_impulse_between_particles_sub(k, i, j+1)
            self.add_impulse_between_particles_sub(k, i+1, j-1)
            self.add_impulse_between_particles_sub(k, i+1, j)
            self.add_impulse_between_particles_sub(k, i+1, j+1)

    def add_impulse_between_bound_particles(self):
        # add impulse between bound pairs
        for k in xrange(self.n):
            for l in self.bonds[k]:
                if k >= l: continue
                rx = self.pos[k,:] - self.pos[l,:]
                rv = self.vel[k,:] - self.vel[l,:]
                d2 = np.sum(rx*rx)
                self.add_impulse_between_bound_pair(k, l, rx, rv, d2)

    def compute_impulse(self):
        self.init_bucket()
        self.add_impulse_from_walls()
        self.add_impulse_between_particles()
        self.add_impulse_between_bound_particles()
        self.pos += self.vel*self.dt

    def update(self):
        self.compute_impulse()
        self.compute_impulse()
        self.compute_impulse()

    def total_energy(self):
        return np.sum(self.mass * np.sum(self.vel**2,axis=1))/2

class HChemViewer:
    INTERVAL = [1000, 100, 20] # in ms
    TEXT_COLOR = (0, 0, 0)
    BACKGROUND_COLOR = (255, 255, 255)
    BOND_COLOR = (0, 0, 0)
    COLORTABLE = [
        (0, 0, 0),
        (255, 0, 0),
        (0, 255, 0),
        (0, 0, 255),
        (255, 0, 255),
        (255, 255, 0),
        (0, 255, 255),
        (192, 192, 192),
        (128, 0, 0),
        (0, 128, 0),
        (0, 0, 128),
        (128, 128, 128),
        (128, 0, 128),
        (128, 128, 0),
        (0, 128, 128),
    ]

    def __init__(self, sim):
        self.sim  = sim
        self.speed = 1
        self.play     = True
        self.stepwise = False
        self.dragged  = None
        self.moving   = False # True when moving a particle
        self.binding  = False # True when binding particles
        self.display_types = False
        self.prev_lclick   = time.time()

        self.root = tk.Tk()
        self.root.title('Tim Hutton\'s Artificial Chemistry Simulator')
        self.root.option_add('*font' 'FixedSys', 14)

        # Infomation variables
        self.show_types = tk.BooleanVar()
        self.show_types.set(True)
        self.show_unbound = tk.BooleanVar()
        self.show_unbound.set(True)
        self.num_particles = tk.IntVar()
        self.num_particles.set(sim.n)
        self.energy = tk.DoubleVar()
        self.energy.set(sim.total_energy())

        self.setup_colormap()
        self.setup_menu()
        self.setup_screen()
        self.setup_statusbar()

    # Setup map table from types to colors
    def setup_colormap(self):
        table = self.COLORTABLE
        self.cmap = {}
        for i in range(len(self.sim.types)):
            self.cmap[self.sim.types[i]] = table[i % len(table)]

    def setup_screen(self):
        self.embed = tk.Frame(self.root,
                width = self.sim.w, height = self.sim.h)
        self.embed.pack()
        os.environ['SDL_WINDOWID'] = str(self.embed.winfo_id())
        os.environ['SDL_VIDEODRIVER'] = 'windib'
        pygame.init()
        self.screen = pygame.display.set_mode((self.sim.w, self.sim.h)
            #, pygame.DOUBLEBUF | pygame.FULLSCREEN | pygame.HWSURFACE
            )
        self.fontsize = 18
        self.font = pygame.font.SysFont(None, self.fontsize)

    def setup_menu(self):
        self.menu = tk.Menu(self.root)
        self.root.configure(m = self.menu)

        file = tk.Menu(self.menu, tearoff = False)
        self.menu.add_cascade(label='File', underline=0, menu=file)
        file.add_command(label='Load Configuration', command=self.load_conf)
        file.add_command(label='Load Rule', command=self.load_rule)
        file.add_command(label='Save Configuration', command=self.save_conf)
        file.add_command(label='Save Rule', command=self.save_rule)
        file.add_command(label='Quit', underline=0,
            command=self.quit, accelerator='Q')
        self.root.bind_all('<Key-q>', self.quit)

        run = tk.Menu(self.menu, tearoff = False)
        self.menu.add_cascade(label='Run', underline=0, menu=run)
        run.add_command(label='Play/Pause', underline=0,
            command=self.cmd_play, accelerator='P')
        self.root.bind_all('<Key-p>', self.cmd_play)
        run.add_command(label='Step execution', underline=1,
            command=self.cmd_stepwise, accelerator='S')
        self.root.bind_all('<Key-s>', self.cmd_stepwise)
        run.add_command(label='Faster', underline=0,
            command=self.cmd_faster, accelerator=u'↑')
        self.root.bind_all('<Key-Up>', self.cmd_faster)
        run.add_command(label='Slower', underline=0,
            command=self.cmd_slower, accelerator=u'↓')
        self.root.bind_all('<Key-Down>', self.cmd_slower)

        show = tk.Menu(self.menu, tearoff = False)
        self.menu.add_cascade(label='Show', underline=0, menu=show)
        show.add_checkbutton(label='Types', variable=self.show_types)
        show.add_checkbutton(label='Unbound Particles',
            variable=self.show_unbound)

    def setup_statusbar(self):
        statusbar = tk.Frame(self.root)
        tk.Label(statusbar, text='# of particles = ').pack(side=tk.LEFT)
        tk.Label(statusbar, textvariable=self.num_particles).pack(side=tk.LEFT)
        tk.Label(statusbar, text=', energy = ').pack(side=tk.LEFT)
        tk.Label(statusbar, textvariable=self.energy).pack(side=tk.LEFT)
        statusbar.pack()

    def cmd_play(self, *event):
        self.play = not self.play

    def cmd_stepwise(self, *event):
        self.play = True
        self.stepwise = True

    def cmd_faster(self, *event):
        self.speed = min(self.speed + 1, len(self.INTERVAL)-1)
    def cmd_slower(self, *event):
        self.speed = max(self.speed - 1, 0)

    def quit(self, *event):
        # XXX: Check whether this simulation is saved or not.
        self.root.destroy()
        sys.exit()

    def load_conf(self, *event):
        print 'Load Conf'

    def load_rule(self, *event):
        print 'Load Rule'

    def save_conf(self, *event):
        print 'Save Conf'

    def save_rule(self, *event):
        print 'Save Rule'

    def check_event(self):
        for event in pygame.event.get():
            if not self.dragged and event.type == pygame.KEYDOWN:
                key = pygame.key.get_pressed()
            elif not self.dragged and event.type == pygame.MOUSEBUTTONDOWN:
                self.play = False
                clicked = self.get_clicked()
                l,m,r = pygame.mouse.get_pressed()

                # Detect double click
                t = time.time()
                double_click = False
                if l and t - self.prev_lclick < 1.0/3:
                    double_click = True
                self.prev_lclick = t

                if double_click:
                    t = inputbox(self.root, 'type')
                    try:
                        self.sim.types[clicked] = self.sim.get_index(t)
                    except:
                        pass
                elif clicked:
                    self.dragged = clicked
                    if l: self.move = True
                    elif r: self.binding = True
            elif self.dragged and event.type == pygame.MOUSEMOTION:
                l,m,r = pygame.mouse.get_pressed()
                if l:
                    self.sim.pos[self.dragged,:] = pygame.mouse.get_pos()
            elif self.dragged and event.type == pygame.MOUSEBUTTONUP:
                if self.binding:
                     clicked = self.get_clicked()
                     if clicked and self.dragged != clicked and\
                             not (clicked in self.sim.bonds[self.dragged]):
                         self.sim.bonds[self.dragged].append(clicked)
                         self.sim.bonds[clicked].append(self.dragged)
                self.move    = False
                self.binding   = False
                self.dragged = None
            elif event.type == pygame.QUIT:
                sys.exit()

    def get_clicked(self):
        for k in xrange(self.sim.n):
            d2 = np.sum((self.sim.pos[k,:] - pygame.mouse.get_pos())**2)
            if d2 < self.sim.r**2:
                return k
                break
        return None

    def loop(self):
        self.draw()
        self.root.mainloop()

    def draw(self):
        sim = self.sim
        if self.play:
            sim.update()
        if self.stepwise:
            self.play = False
            self.stepwise = False

        # Update information
        self.energy.set(sim.total_energy())

        self.screen.fill(self.BACKGROUND_COLOR)

        # Draw particles
        for k in xrange(sim.n):
            if not self.show_unbound.get() and not sim.bonds[k]:
                continue
            pygame.draw.circle(self.screen, self.cmap[sim.typ[k]],
                (int(sim.pos[k,0]), int(sim.pos[k,1])), sim.r, 1)
        if self.show_types.get():
            for k in xrange(sim.n):
                if not self.show_unbound.get() and not sim.bonds[k]:
                    continue
                t = sim.typ[k] + str(sim.state[k])
                text = self.font.render(t, False, self.TEXT_COLOR)
                rect = text.get_rect()
                rect.centerx = sim.pos[k, 0]
                rect.centery = sim.pos[k, 1]
                self.screen.blit(text, rect)

        # Draw bonds
        for k in xrange(sim.n):
            for l in sim.bonds[k]:
                if k >= l: continue
                pygame.draw.line(self.screen, self.BOND_COLOR,
                    self.sim.pos[k,:], self.sim.pos[l,:])

        pygame.display.flip()
        self.root.after(self.INTERVAL[self.speed], self.draw)

#            # Other info
#            if self.binding:
#                pygame.draw.line(self.screen, self.BOND_COLOR,
#                    self.sim.pos[self.dragged,:],
#                    pygame.mouse.get_pos())
#
#            self.check_event()
#            pygame.display.update()
#            self.root.update()

#class HChemViewer:
#    INFO = [
#    "(P) play/pause, (F) stepwise, (Q) quit, (S) save snapshot, (L) load snapshot, (R) load rule, (T) show/hide particle types, (up) speed up, (down) slow down",
#    "(left drag) move particle, (right drag) bind particles, (double click) change particle type'
#    ]
#
#    def __init__(self, sim, interval=0.1):
#        self.sim      = sim
#        self.w        = sim.w
#        self.h        = sim.h
#        self.interval = interval
#
#        self.pause    = False
#        self.stepwise = False
#        self.dragged  = None
#        self.moving   = False # True when moving a particle
#        self.binding  = False # True when binding particles
#        self.display_types = False
#        self.prev_lclick   = time.time()
#
#        self.setup_colormap()
#        self.setup_screen()
#        self.setup_info_text()
#
#    # Setup map table from types to colors
#    def setup_colormap(self):
#        table = self.COLORTABLE
#        self.cmap = {}
#        for i in range(len(self.sim.types)):
#            self.cmap[self.sim.types[i]] = table[i % len(table)]
#
#    def setup_screen(self):
#        pygame.init()
#        self.screen = pygame.display.set_mode((self.w, self.h)
#            #, pygame.DOUBLEBUF | pygame.FULLSCREEN | pygame.HWSURFACE
#            )
#        pygame.display.set_caption('Tim Hutton's Artificial Chemistry')
#        self.fontsize = 18
#        self.font = pygame.font.SysFont(None, self.fontsize)
#
#    def setup_info_text(self):
#        info_texts = self.INFO + map(str, self.sim.rules)
#        self.info = map(lambda text:
#                self.font.render(text, False, self.TEXT_COLOR),
#                info_texts)
