// Copyright 2012 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/**
 * Island represents a grouping of objects to apply DFS solving
 * for movement and collisions.
 */

part of box2d;

class Island {
  ContactListener listener;

  List<Body> bodies = <Body>[];
  List<Contact> contacts = <Contact>[];
  List<Joint> joints = <Joint>[];

  int bodyCount = 0;
  int jointCount = 0;
  int contactCount = 0;

  int positionIterationCount;

  // Pool objects to cut down on object creation.
  ContactSolver _contactSolver = new ContactSolver();
  Vector2 _translation = new Vector2.zero();
  ContactImpulse impulse = new ContactImpulse();

  Island();

  //TODO(gregbglw): No need to keep capacity, count and array for these items as
  // in C. Simply measure the length of the array, for example, in order to
  // determine capacity.
  void init(int argBodyCapacity, int argContactCapacity, int argJointCapacity,
      ContactListener argListener) {
    bodyCount = 0;
    contactCount = 0;
    jointCount = 0;

    listener = argListener;

    bodies.length = argBodyCapacity;
    contacts.length = argContactCapacity;
    joints.length = argJointCapacity;
  }

  void clear() {
    bodyCount = 0;
    contactCount = 0;
    jointCount = 0;
  }

  void solve(TimeStep step, Vector2 gravity, bool allowSleep) {
    // Integrate velocities and apply damping.
    for (int i = 0; i < bodyCount; ++i) {
      Body b = bodies[i];

      if (b.type != BodyType.DYNAMIC) {
        continue;
      }

      final velocityDelta = new Vector2(
          (b._force.x * b.invMass + gravity.x) * step.dt,
          (b._force.y * b.invMass + gravity.y) * step.dt);
      b.linearVelocity.add(velocityDelta);
      num newAngularVelocity = b.angularVelocity +
          (step.dt * b.invInertia * b._torque);
      b.angularVelocity = newAngularVelocity;

      num a = (1.0 - step.dt * b.linearDamping);
      num a1 = (0.0 > (a < 1.0 ? a : 1.0) ? 0.0 : (a < 1.0 ? a : 1.0));
      b.linearVelocity.scale(a1);

      num a2 = (1.0 - step.dt * b.angularDamping);
      num b1 = (a2 < 1.0 ? a2 : 1.0);
      b.angularVelocity *= 0.0 > b1 ? 0.0 : b1;
    }

    // Partition contacts so that contacts with static bodies are solved last.
    int i1 = -1;
    for (int i2 = 0; i2 < contactCount; ++i2) {
      Fixture fixtureA = contacts[i2].fixtureA;
      Fixture fixtureB = contacts[i2].fixtureB;
      Body bodyA = fixtureA.body;
      Body bodyB = fixtureB.body;
      bool nonStatic = bodyA.type != BodyType.STATIC && bodyB.type
          != BodyType.STATIC;
      if (nonStatic){
        ++i1;
        //Swap(contacts[i1], contacts[i2]);
        Contact temp = contacts[i1];
        contacts[i1] = contacts[i2];
        contacts[i2] = temp;
      }
    }

    // Initialize velocity constraints.
    _contactSolver.init(contacts, contactCount, step.dtRatio);
    _contactSolver.warmStart();

    for (int i = 0; i < jointCount; ++i) {
      joints[i].initVelocityConstraints(step);
    }

    for (int i = 0; i < step.velocityIterations; ++i) {
      for (int j = 0; j < jointCount; ++j) {
        joints[j].solveVelocityConstraints(step);
      }
      _contactSolver.solveVelocityConstraints();
    }

    // Post-solve (store impulses for warm starting).
    _contactSolver.storeImpulses();

    // Integrate positions.
    Vector2 temp = new Vector2.zero();
    for (int i = 0; i < bodyCount; ++i) {
      Body b = bodies[i];

      if (b.type == BodyType.STATIC) {
        continue;
      }

      // Check for large velocities.
      _translation.setFrom(b.linearVelocity);
      _translation.scale(step.dt);
      if (_translation.dot(_translation) > Settings.MAX_TRANSLATION_SQUARED) {
        num ratio = Settings.MAX_TRANSLATION / _translation.length;
        b.linearVelocity.scale(ratio);
      }

      num rotation = step.dt * b.angularVelocity;
      if (rotation * rotation > Settings.MAX_ROTATION_SQUARED) {
        num ratio = Settings.MAX_ROTATION / rotation.abs();
        b.angularVelocity *= ratio;
      }

      // Store positions for continuous collision.
      b.sweep.centerZero.setFrom(b.sweep.center);
      b.sweep.angleZero = b.sweep.angle;

      // Integrate
      temp.setFrom(b.linearVelocity);
      temp.scale(step.dt);
      b.sweep.center.add(temp);
      b.sweep.angle += step.dt * b.angularVelocity;

      // Compute new transform
      b.synchronizeTransform();

      // Note: shapes are synchronized later.
    }

    // Iterate over constraints.
    for (int i = 0; i < step.positionIterations; ++i){
      bool contactsOkay =
          _contactSolver.solvePositionConstraints(Settings.CONTACT_BAUMGARTE);

      bool jointsOkay = true;
      for (int j = 0; j < jointCount; ++j){
        bool jointOkay =
            joints[j].solvePositionConstraints(Settings.CONTACT_BAUMGARTE);
        jointsOkay = jointsOkay && jointOkay;
      }

      if (contactsOkay && jointsOkay) {
        // Exit early if the position errors are small.
        break;
      }
    }

    report(_contactSolver.constraints);

    _contactSolver.clearReferences();


    if (allowSleep) {
      num minSleepTime = Settings.BIG_NUMBER;

      num linTolSqr = Settings.LINEAR_SLEEP_TOLERANCE
          * Settings.LINEAR_SLEEP_TOLERANCE;
      num angTolSqr = Settings.ANGULAR_SLEEP_TOLERANCE
          * Settings.ANGULAR_SLEEP_TOLERANCE;

      for (int i = 0; i < bodyCount; ++i) {
        Body b = bodies[i];
        if (b.type == BodyType.STATIC) {
          continue;
        }

        if ((b.flags & Body.AUTO_SLEEP_FLAG) == 0) {
          b.sleepTime = 0.0;
          minSleepTime = 0.0;
        }

        if ((b.flags & Body.AUTO_SLEEP_FLAG) == 0 ||
            b.angularVelocity * b.angularVelocity > angTolSqr ||
            b.linearVelocity.dot(b.linearVelocity) > linTolSqr){
          b.sleepTime = 0.0;
          minSleepTime = 0.0;
        } else {
          b.sleepTime += step.dt;
          minSleepTime = Math.min(minSleepTime, b.sleepTime);
        }
      }

      if (minSleepTime >= Settings.TIME_TO_SLEEP) {
        for (int i = 0; i < bodyCount; ++i) {
          Body b = bodies[i];
          b.awake = false;
        }
      }
    }

  }

  /** Adds a body to the Island. */
  void addBody(Body body){
    assert(bodyCount < bodies.length);
    body.islandIndex = bodyCount;
    bodies[bodyCount++] = body;
  }

  /** Add a contact to the Island. */
  void addContact(Contact contact) {
    assert(contactCount < contacts.length);
    contacts[contactCount++] = contact;
  }

  /** Add a joint to the Island. */
  void addJoint(Joint joint) {
    assert(jointCount < joints.length);
    joints[jointCount++] = joint;
  }

  void report(List<ContactConstraint> constraints) {
    if (listener == null) {
      return;
    }

    for (int i = 0; i < contactCount; ++i) {
      Contact c = contacts[i];

      ContactConstraint cc = constraints[i];

      for (int j = 0; j < cc.pointCount; ++j) {
        impulse.normalImpulses[j] = cc.points[j].normalImpulse;
        impulse.tangentImpulses[j] = cc.points[j].tangentImpulse;
      }

      listener.postSolve(c, impulse);
    }
  }
}

/**
 * This is an internal structure
 */
class Position {
  Vector2 x;
  num a;

  Position() {
    x = new Vector2.zero();
    a = 0;
  }
}

/**
 * This is an internal structure
 */
class Velocity {
  Vector2 v;
  num a;

  Velocity() {
    v = new Vector2.zero();
    a = 0;
  }
}
