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
 *  Used to warm start distance.
 */

part of box2d;

class SimplexCache {
  /** length or area */
  double metric = 0.0;

  int count = 0;

  /** vertices on shape A */
  final List<int> indexA =
    new List<int>.generate(3, (i) => Settings.MAX_INTEGER);

  /** vertices on shape B */
  final List<int> indexB =
    new List<int>.generate(3, (i) => Settings.MAX_INTEGER);

  /**
   * Constructs a new SimplexCache.
   */
  SimplexCache();

  /**
   * Sets this cache equal to the given cache.
   */
  void setFrom(SimplexCache sc) {
    indexA.setRange(0, indexA.length, sc.indexA);
    indexB.setRange(0, indexB.length, sc.indexB);
    metric = sc.metric;
    count = sc.count;
  }
}
